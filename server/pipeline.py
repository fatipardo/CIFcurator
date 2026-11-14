#!/usr/bin/env python
# coding: utf-8
#
r"""Run the whole CIFcurator workflow (Steps 1 to 5) without asking any question.

This is what the CIFcurator server runs for every job, but it can also be run on its own:

  python server/pipeline.py "triosephosphate isomerase" --projdir ./TIM

Steps, all written inside projdir:
  1  CIFcurator_Download.py: download the structures matching the keyword (raw_bank) and clean
     them (clean_bank). With --csm, the computed structure models (CSMs) of AlphaFold DB that
     match the keyword are downloaded too (at most --max-models of them); a CSM has ONE PROTEIN
     CHAIN ONLY (no nucleic acids, ligands, waters or other copies).
     With --exclude-file, the structures it names are left out of Steps 3-5 (they stay in
     clean_bank); the structures left out are listed in projdir/excluded_structures.txt.
  3  MolID: cluster the chains with MMseqs2 (clustering_storage/), then name the chains with
     --input-file if given, otherwise with the input file generated from the clusters
     (input_file.txt, menu option 4 of CIFcurator_MolID.py). Concatenations are all accepted.
     Output: standard_MolID_bank.
  4  CIFcurator_ChainStandardization.py: Standard Sequences from --reference-json if given,
     otherwise generated from all the structures (consensus). Output: standard_ChainID_bank
     (with StandardSequences.json and ChainStandardizationRecord.txt).
  5  CIFcurator_ResidueStandardization.py: residue renumbering with the occupancy threshold
     (--threshold, default 80). Output: standard_ResidueID_bank (with OldResID_NewResID_Map.csv).
     The MUSCLE alignments are kept in residue_alignments/.

The output of every step goes to projdir/logs/. The PDB sequence file (pdb_seqres.fasta, about
400 MB) is NOT written to projdir: it is kept in --seqfile and shared between runs.

Exit code 0 when every step worked. Otherwise the message of the failed step is printed and
written to projdir/logs/error.txt.

The steps are run as separate processes, with projdir (Step 3) or projdir/residue_alignments
(Step 5) as their working directory, because Step 3 changes the working directory while it
clusters and Step 5 writes its alignment files in the working directory.
"""

import argparse, glob, json, os, re, subprocess, sys, tempfile, time

SERVER_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_DIR = os.path.dirname(SERVER_DIR)
SCRIPTS_DIR = os.path.join(REPO_DIR, 'scripts')
SRC_DIR = os.path.join(REPO_DIR, 'src')

DEFAULT_THRESHOLD = 80
# The PDB releases new entries, and a new pdb_seqres file, every Wednesday at 00:00 UTC
# (https://www.rcsb.org/docs/general-help/data-release-schedule). The sequence file is
# downloaded again when it is older than the last release.
RELEASE_WEEKDAY = 2  # Wednesday (Monday is 0)

# Name of every step, in the order they run. The server shows them as they are written to
# projdir/logs/steps.json.
STEPS = [
    ('download', 'Step 1: download and clean structures'),
    ('exclude', 'Leave out the structures you listed'),
    ('molid', 'Step 3: entity (MolID) chain naming'),
    ('chainid', 'Step 4: chain ID standardization'),
    ('residueid', 'Step 5: residue number standardization'),
    ('zip', 'Package results'),
]


class StepFailed(Exception):
    pass


def package_env():
    """
    Environment in which `import CIFcurator` gives this checkout's src/, even if another
    version of CIFcurator is installed (the package lives in src/, see setup.py package_dir).
    """
    shim = os.path.join(tempfile.gettempdir(), 'cifcurator_server_pylib_{0}'.format(os.getuid()))
    link = os.path.join(shim, 'CIFcurator')
    os.makedirs(shim, exist_ok=True)
    if os.path.realpath(link) != os.path.realpath(SRC_DIR):
        tmp = link + '.{0}'.format(os.getpid())
        os.symlink(SRC_DIR, tmp)
        os.replace(tmp, link)
    env = os.environ.copy()
    env['PYTHONPATH'] = os.pathsep.join(p for p in (shim, env.get('PYTHONPATH')) if p)
    env['PYTHONWARNINGS'] = 'ignore::SyntaxWarning'
    env['PYTHONUNBUFFERED'] = '1'
    return env


class Pipeline:
    def __init__(self, keyword, projdir, seqfile, input_file=None, reference_json=None,
                 threshold=DEFAULT_THRESHOLD, csm=False, max_models=None, exclude=None,
                 update_seqfile=False, seqres_max_age_days=None):
        self.keyword = keyword
        self.exclude = list(exclude or [])  # names given by the user (see match_structures)
        self.excluded_files = []  # the files of clean_bank they matched
        self.csm = csm
        self.max_models = max_models
        self.projdir = os.path.abspath(projdir)
        self.seqfile = os.path.abspath(seqfile)
        self.input_file = os.path.abspath(input_file) if input_file else None
        self.reference_json = os.path.abspath(reference_json) if reference_json else None
        self.threshold = threshold
        self.update_seqfile = update_seqfile
        self.seqres_max_age_days = seqres_max_age_days
        self.logdir = os.path.join(self.projdir, 'logs')
        self.env = package_env()
        self.status = {key: {'title': title, 'state': 'pending'} for key, title in STEPS}

    # -- bookkeeping ------------------------------------------------------------------
    def write_status(self):
        path = os.path.join(self.logdir, 'steps.json')
        with open(path + '.tmp', 'w') as f:
            json.dump({'steps': [dict(key=k, **self.status[k]) for k, _ in STEPS]}, f, indent=2)
        os.replace(path + '.tmp', path)

    def set_state(self, key, state, **extra):
        self.status[key]['state'] = state
        self.status[key].update(extra)
        self.write_status()

    def run(self, key, cmd, cwd, answers=None):
        """Run one step as a process, its output going to logs/<key>.log."""
        os.makedirs(cwd, exist_ok=True)
        log_path = os.path.join(self.logdir, key + '.log')
        with open(log_path, 'w') as log:
            log.write('$ ' + ' '.join(cmd) + '\n')
            if answers is not None:
                log.write('# answers on standard input: {0!r}\n'.format(answers))
            log.flush()
            proc = subprocess.Popen(cmd, cwd=cwd, env=self.env, stdout=log, stderr=subprocess.STDOUT,
                                    stdin=subprocess.PIPE if answers is not None else subprocess.DEVNULL, text=True)
            proc.communicate(answers)
        if proc.returncode != 0:
            raise StepFailed('{0} failed (exit code {1}); see logs/{2}.log:\n{3}'.format(
                self.status[key]['title'], proc.returncode, key, tail(log_path)))

    def step(self, key, function):
        self.set_state(key, 'running', started=time.time())
        try:
            function()
        except Exception as err:
            self.set_state(key, 'failed', finished=time.time(), error=str(err))
            raise
        self.set_state(key, 'done', finished=time.time())

    # -- the steps ----------------------------------------------------------------------
    def download(self):
        cmd = [sys.executable, os.path.join(SCRIPTS_DIR, 'CIFcurator_Download.py'), self.keyword,
               '--projdir', self.projdir, '--seqfile', self.seqfile]
        if self.csm:
            cmd += ['--csm', 'True']
            if self.max_models is not None:
                cmd += ['--max-models', str(self.max_models)]
        if self.update_seqfile or seqfile_is_stale(self.seqfile, self.seqres_max_age_days):
            cmd.append('--update')
        self.run('download', cmd, cwd=self.projdir)
        if not glob.glob(os.path.join(self.projdir, 'clean_bank', '*.cif')):
            raise StepFailed('No structure was downloaded for the keyword {0!r}.'.format(self.keyword))

    def exclude_structures(self):
        """Find the files of clean_bank named in self.exclude; they are not given to Step 3."""
        log_path = os.path.join(self.logdir, 'exclude.log')
        files = [os.path.basename(f) for f in glob.glob(os.path.join(self.projdir, 'clean_bank', '*.cif'))]
        matched, not_found = match_structures(self.exclude, files)
        self.excluded_files = matched
        kept = len(files) - len(matched)
        lines = ['Structures left out of Steps 3-5 (they are still in clean_bank):']
        lines += ['  ' + f for f in matched] or ['  none']
        if not_found:
            lines += ['', 'WARNING: these names match no downloaded structure, so nothing was left out for them:']
            lines += ['  ' + name for name in not_found]
        lines += ['', '{0} of the {1} structures of clean_bank go on to Step 3.'.format(kept, len(files))]
        text = '\n'.join(lines) + '\n'
        for path in (log_path, os.path.join(self.projdir, 'excluded_structures.txt')):
            with open(path, 'w') as f:
                f.write(text)
        print(text)
        if kept == 0:
            raise StepFailed('Every downloaded structure was left out: Steps 3-5 need at least one.')
        summary = '{0} left out, {1} kept'.format(len(matched), kept)
        if not_found:
            summary += '; not found: ' + ', '.join(not_found[:10]) + (' ...' if len(not_found) > 10 else '')
        self.status['exclude']['summary'] = summary

    # The Step 3-5 scripts write into target_dir but do not create it
    def molid(self):
        os.makedirs(os.path.join(self.projdir, 'standard_MolID_bank'), exist_ok=True)
        cmd = [sys.executable, os.path.abspath(__file__), '_molid',
               os.path.join(self.projdir, 'clean_bank'), os.path.join(self.projdir, 'standard_MolID_bank'),
               self.projdir, self.seqfile]
        if self.input_file:
            cmd += ['--input-file', self.input_file]
        if self.excluded_files:
            path = os.path.join(self.logdir, 'excluded_files.txt')
            with open(path, 'w') as f:
                f.write('\n'.join(self.excluded_files) + '\n')
            cmd += ['--exclude-file', path]
        self.run('molid', cmd, cwd=self.projdir)

    def chainid(self):
        os.makedirs(os.path.join(self.projdir, 'standard_ChainID_bank'), exist_ok=True)
        if self.reference_json:
            answers = '2\n5\n{0}\n4\n4\n'.format(self.reference_json)
        else:
            answers = '2\n4\n4\n4\n'
        cmd = [sys.executable, os.path.join(SCRIPTS_DIR, 'CIFcurator_ChainStandardization.py'),
               os.path.join(self.projdir, 'standard_MolID_bank'), os.path.join(self.projdir, 'standard_ChainID_bank')]
        self.run('chainid', cmd, cwd=self.projdir, answers=answers)

    def residueid(self):
        os.makedirs(os.path.join(self.projdir, 'standard_ResidueID_bank'), exist_ok=True)
        answers = '1\n4\n{0}\n4\n3\nQUIT\n'.format(self.threshold)
        cmd = [sys.executable, os.path.join(SCRIPTS_DIR, 'CIFcurator_ResidueStandardization.py'),
               os.path.join(self.projdir, 'standard_ChainID_bank'), os.path.join(self.projdir, 'standard_ResidueID_bank')]
        self.run('residueid', cmd, cwd=os.path.join(self.projdir, 'residue_alignments'), answers=answers)
        if not glob.glob(os.path.join(self.projdir, 'standard_ResidueID_bank', '*.cif')):
            raise StepFailed('Step 5 wrote no CIF file; see logs/residueid.log:\n'
                             + tail(os.path.join(self.logdir, 'residueid.log')))

    def run_all(self, zip_path=None):
        os.makedirs(self.logdir, exist_ok=True)
        with open(os.path.join(self.logdir, 'parameters.json'), 'w') as f:
            json.dump({'keyword': self.keyword, 'threshold': self.threshold,
                       'csm': self.csm, 'max_models': self.max_models, 'exclude': self.exclude,
                       'input_file': self.input_file, 'reference_json': self.reference_json,
                       'seqfile': self.seqfile, 'cifcurator_checkout': REPO_DIR,
                       'cifcurator_commit': git_commit()}, f, indent=2)
        self.write_status()
        try:
            self.step('download', self.download)
            if self.exclude:
                self.step('exclude', self.exclude_structures)
            else:
                self.set_state('exclude', 'skipped')
            self.step('molid', self.molid)
            self.step('chainid', self.chainid)
            self.step('residueid', self.residueid)
        except Exception as err:
            with open(os.path.join(self.logdir, 'error.txt'), 'w') as f:
                f.write(str(err) + '\n')
            # The files written so far are still worth having
            if zip_path:
                self.step('zip', lambda: make_zip(self.projdir, zip_path))
            raise
        if zip_path:
            self.step('zip', lambda: make_zip(self.projdir, zip_path))
        else:
            self.set_state('zip', 'skipped')


# -- Step 3 without menus ---------------------------------------------------------------
def run_molid(source_dir, target_dir, projdir, seqfile, input_file=None, exclude=()):
    """
    Step 3, as CIFcurator_MolID.py does it with the answers 4 (generate input file),
    3 (enter it), 7 (continue) and, if there are concatenations, 5 (accept all), then 6
    (finalize). With input_file, that file is entered instead of the generated one
    (clustering still runs, as in the script, and clustering_storage/ is kept).
    The files of source_dir named in exclude (file names) are left out.
    """
    from CIFcurator import cifcuratormolidcifutils as molidutils
    from CIFcurator import cifcluster

    exclude = set(exclude)
    filelist = sorted(f for f in glob.glob(source_dir + '/*.cif') if os.path.basename(f) not in exclude)
    if exclude:
        print('Leaving out {0} structures: {1}'.format(len(exclude), ', '.join(sorted(exclude))))
    if not filelist:
        sys.exit('No CIF file in {0}'.format(source_dir) + (' after leaving out the excluded structures' if exclude else ''))
    cifcluster.run_clustering_pipeline(filelist, projdir, seqfile)
    cluster_pdbchids, outlier_case_entries = cifcluster.finding_problematic_cluster_results(projdir + '/clustering_storage/')

    master_molID_class_list = molidutils.pdb_to_masterlist(filelist)
    unique_molID_occur_map = molidutils.CreateMasterUniqueMolIDMap(master_molID_class_list)
    molIDConversion_list = molidutils.uniquelist_to_conversionlist(unique_molID_occur_map)
    MolID_to_files_map = molidutils.CreateMasterUniqueMolIDMapWithFileName(master_molID_class_list)
    MolID_occur_dict_of_lists = molidutils.CreateMasterUniqueMolIDOccursLIST(master_molID_class_list)
    MolID_ChainID_dict_of_lists = molidutils.CreateMasterUniqueMolIDinitialChainIDsLIST(master_molID_class_list)

    if input_file is None:
        print('Generating the input file from the MMseqs2 clusters (menu option 4)')
        cifcluster.generate_chain_assignments(molIDConversion_list, master_molID_class_list, cluster_pdbchids,
                                              outlier_case_entries, projdir)
        input_file = os.path.join(projdir, 'input_file.txt')
    print('Entering the input file {0} (menu option 3)'.format(input_file))
    user_molID_chID_map = molidutils.read_input_file(input_file)
    for molIDConversion in molIDConversion_list:  # as molidutils.add_user_conversion
        for key in user_molID_chID_map:
            if str(key) == str(molIDConversion.molID):
                molIDConversion.add_chID_list(user_molID_chID_map[key])
        molIDConversion.check_for_completeness()

    if molidutils.check_complete(molIDConversion_list) != '1':
        print('\nThese entities are missing chain IDs:')
        molidutils.show_unassigned_conversion(molIDConversion_list)
        sys.exit('ERROR: the input file {0} does not give every entity all its chain IDs '
                 '(the list above shows what is missing; the entity names must match exactly, '
                 'and there must be one chain ID per copy).'.format(os.path.basename(input_file)))

    print('Saving track changes (menu option 7)')
    molidutils.show_full_conversion_and_file_list_by_number_chains(MolID_ChainID_dict_of_lists, molIDConversion_list,
                                                                   MolID_to_files_map, MolID_occur_dict_of_lists, target_dir)
    master_molID_class_list = molidutils.update_masterlist(master_molID_class_list, molIDConversion_list)
    if molidutils.problem_counter(master_molID_class_list) != 0:
        print('Accepting ALL concatenations (concatenation menu option 5)')
        for newchain in molidutils.return_unassigned_conversion(master_molID_class_list, step='concatenation'):
            master_molID_class_list = molidutils.list_accept_concatenations_auto(master_molID_class_list, newchain)[0]
        left = molidutils.problem_counter(master_molID_class_list)
        if left != 0:
            sys.exit('ERROR: {0} concatenations could not be accepted'.format(left))
    print('Finalizing curation (menu option 6)')
    molidutils.masterlist_to_pdb(filelist, master_molID_class_list, target_dir=target_dir)
    print('Done')


# -- helpers --------------------------------------------------------------------------------
def read_name_list(text):
    """Names separated by new lines, commas, semicolons or spaces; '#' starts a comment."""
    names = []
    for line in text.splitlines():
        for name in re.split(r'[\s,;]+', line.split('#', 1)[0]):
            if name and name not in names:
                names.append(name)
    return names


def structure_keys(name):
    """
    The ways a structure can be named, lowercase: the file name with and without .cif and, for
    the files of an entry (pdb_00001abc+01.cif, AF-P00780-F1+01.cif), the entry without the
    assembly number (pdb_00001abc, af-p00780-f1) and the 4-character PDB ID (1abc).
    """
    name = os.path.basename(name.strip()).lower()
    stem = name[:-4] if name.endswith('.cif') else name
    keys = {name, stem}
    entry = stem.split('+', 1)[0]
    keys.add(entry)
    match = re.match(r'^pdb_0000([0-9a-z]{4})$', entry)
    if match:
        keys.add(match.group(1))
    return keys


def match_structures(names, files):
    """
    The files (names of the files of clean_bank) matched by the names the user gave, and the names
    that match no file. A name matches a file when it is its file name (with or without .cif), its
    entry (all the assemblies of an entry: pdb_00001abc, AF-P00780-F1), or its PDB ID (1abc, or
    1ABC: case does not matter).
    """
    matched, not_found = [], []
    for name in names:
        query = name.strip().lower()
        query = query[:-4] if query.endswith('.cif') else query
        if re.match(r'^[0-9][0-9a-z]{3}$', query):  # a PDB ID: its assemblies are pdb_0000<id>+NN.cif
            query = 'pdb_0000' + query
        hits = [f for f in files if query in structure_keys(f)]
        if not hits:
            not_found.append(name)
        matched += [f for f in hits if f not in matched]
    return sorted(matched), not_found


def last_pdb_release(now=None):
    """Time (seconds since the epoch) of the last weekly PDB release, Wednesday 00:00 UTC."""
    import datetime
    now = datetime.datetime.fromtimestamp(time.time() if now is None else now, datetime.timezone.utc)
    midnight = now.replace(hour=0, minute=0, second=0, microsecond=0)
    return (midnight - datetime.timedelta(days=(now.weekday() - RELEASE_WEEKDAY) % 7)).timestamp()


def seqfile_is_stale(seqfile, max_age_days=None):
    """
    True when the sequence file is missing or empty, or was downloaded before the last weekly PDB
    release (or, with max_age_days, is older than that many days).
    """
    if not os.path.isfile(seqfile) or os.path.getsize(seqfile) == 0:
        return True
    if max_age_days is not None:
        return (time.time() - os.path.getmtime(seqfile)) / 86400 > max_age_days
    return os.path.getmtime(seqfile) < last_pdb_release()


def tail(path, lines=25):
    try:
        with open(path, errors='replace') as f:
            return ''.join(f.readlines()[-lines:])
    except OSError:
        return ''


def git_commit():
    try:
        return subprocess.run(['git', '-C', REPO_DIR, 'rev-parse', '--short', 'HEAD'], capture_output=True,
                              text=True, timeout=10).stdout.strip() or None
    except (OSError, subprocess.SubprocessError):
        return None


def make_zip(projdir, zip_path):
    """Every file of projdir, in a folder named after the zip file."""
    import zipfile
    top = os.path.splitext(os.path.basename(zip_path))[0]
    partial = zip_path + '.part'
    with zipfile.ZipFile(partial, 'w', zipfile.ZIP_DEFLATED) as z:
        for root, dirs, files in os.walk(projdir):
            dirs.sort()
            for name in sorted(files):
                path = os.path.join(root, name)
                if os.path.islink(path) or os.path.abspath(path) == os.path.abspath(partial):
                    continue
                z.write(path, os.path.join(top, os.path.relpath(path, projdir)))
    os.replace(partial, zip_path)


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if argv and argv[0] == '_molid':
        p = argparse.ArgumentParser(prog='pipeline.py _molid')
        p.add_argument('source_dir'); p.add_argument('target_dir'); p.add_argument('projdir'); p.add_argument('seqfile')
        p.add_argument('--input-file')
        p.add_argument('--exclude-file')
        a = p.parse_args(argv[1:])
        exclude = []
        if a.exclude_file:
            with open(a.exclude_file) as f:
                exclude = read_name_list(f.read())
        run_molid(a.source_dir, a.target_dir, a.projdir, a.seqfile, a.input_file, exclude)
        return 0

    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0],
                                     epilog=__doc__.split('\n\n', 1)[1],
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('keyword', help='Keyword to search for in the PDB, e.g. "triosephosphate isomerase" (in quotes)')
    parser.add_argument('--projdir', required=True, help='Project directory (created if needed)')
    parser.add_argument('--seqfile', default=None,
                        help='PDB sequence file, downloaded if missing or older than the last weekly PDB release '
                             '(default: server_data/pdb_seqres.fasta next to this script)')
    parser.add_argument('--seqres-max-age-days', type=float, default=None, metavar='DAYS',
                        help='Download the sequence file again when it is older than DAYS days (default: when it '
                             'was downloaded before the last weekly PDB release, Wednesday 00:00 UTC)')
    parser.add_argument('--update', action='store_true', help='Download the sequence file again now')
    parser.add_argument('--input-file', default=None,
                        help='Step 3 input file (ENTITY NAME:chainID,chainID,... per line). '
                             'Default: generated from the MMseqs2 clusters')
    parser.add_argument('--reference-json', default=None,
                        help='Step 4 Standard Sequences (StandardSequences.json). Default: consensus of all structures')
    parser.add_argument('--threshold', type=int, default=DEFAULT_THRESHOLD,
                        help='Step 5 occupancy threshold for residue renumbering, 1-100 (default: %(default)s)')
    parser.add_argument('--csm', action='store_true',
                        help='Step 1: also download the computed structure models (CSMs) of AlphaFold DB that match '
                             'the keyword. A CSM has one protein chain only')
    parser.add_argument('--max-models', type=int, default=None, metavar='N',
                        help='With --csm, download at most N CSMs, best matches first (default: all of them)')
    parser.add_argument('--exclude-file', default=None, metavar='FILE',
                        help='Leave the structures named in FILE out of Steps 3-5 (they stay in clean_bank). One or '
                             'more names per line (separated by spaces or commas, # starts a comment): a file name '
                             '(pdb_00001abc+01.cif, with or without .cif), an entry (pdb_00001abc or AF-P00780-F1: all '
                             'its assemblies) or a PDB ID (1abc). Default: all the structures go on to Step 3')
    parser.add_argument('--zip', default=None, help='Also write every file of projdir to this zip file')
    a = parser.parse_args(argv)
    if not 1 <= a.threshold <= 100:
        parser.error('--threshold must be between 1 and 100')
    if a.max_models is not None and a.max_models < 1:
        parser.error('--max-models must be at least 1')
    exclude = []
    if a.exclude_file:
        with open(a.exclude_file) as f:
            exclude = read_name_list(f.read())
    seqfile = a.seqfile or os.path.join(SERVER_DIR, 'server_data', 'pdb_seqres.fasta')
    os.makedirs(os.path.dirname(os.path.abspath(seqfile)), exist_ok=True)
    pipeline = Pipeline(a.keyword, a.projdir, seqfile, input_file=a.input_file, reference_json=a.reference_json,
                        threshold=a.threshold, csm=a.csm, max_models=a.max_models, exclude=exclude, update_seqfile=a.update, seqres_max_age_days=a.seqres_max_age_days)
    try:
        pipeline.run_all(zip_path=a.zip)
    except Exception as err:
        print(err, file=sys.stderr)
        return 1
    print('Done. Results in {0}'.format(pipeline.projdir) + (' and {0}'.format(a.zip) if a.zip else ''))
    return 0


if __name__ == '__main__':
    sys.exit(main())
