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

Rerun from a later step (--start-step): the steps before it are not run again, their output
already in projdir is used (state "reused"). With --reuse-from, that output is first taken from the
project directory of an earlier run by hard links (no copy, no extra space; the earlier run is not
changed, as no step writes into the files it reads). Everything else must match the earlier run
(same keyword and computed-model options); the list of structures to leave out is that of the
earlier run unless the rerun starts from the exclusion step.

The output of every step goes to projdir/logs/. The PDB sequence file (pdb_seqres.fasta, about
400 MB) is NOT written to projdir: it is kept in --seqfile and shared between runs.

Exit code 0 when every step worked. Otherwise the message of the failed step is printed and
written to projdir/logs/error.txt.

Every path given to a step is relative to its working directory, so that the files written in
projdir (logs, alignment headers, track-changes records) do not show where projdir is: on the server,
that path holds the job number, and the zip file may be shared.

The steps are run as separate processes, with projdir (Step 3) or projdir/residue_alignments
(Step 5) as their working directory, because Step 3 changes the working directory while it
clusters and Step 5 writes its alignment files in the working directory.
"""

import argparse, glob, json, os, re, shutil, subprocess, sys, tempfile, time

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

# Steps that a rerun can start from, and what every step writes in projdir (what a rerun reuses)
RERUN_STEPS = ['exclude', 'molid', 'chainid', 'residueid']
STEP_OUTPUTS = {
    'download': ['raw_bank', 'clean_bank', 'logs/download.log'],
    'exclude': ['excluded_structures.txt', 'logs/exclude.log', 'logs/excluded_files.txt'],
    'molid': ['clustering_storage', 'input_file.txt', 'standard_MolID_bank', 'logs/molid.log'],
    'chainid': ['standard_ChainID_bank', 'logs/chainid.log'],
    'residueid': ['residue_alignments', 'standard_ResidueID_bank', 'logs/residueid.log'],
}


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
                 update_seqfile=False, seqres_max_age_days=None, start_step=None, max_entries=None):
        self.keyword = keyword
        self.max_entries = max_entries  # None: no limit on the number of PDB entries to download
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
        self.start_step = start_step  # None: from Step 1; else one of RERUN_STEPS
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

    def rel(self, path, cwd=None):
        """path relative to cwd (default: projdir), as the steps are given it."""
        return os.path.relpath(path, cwd or self.projdir)

    def run(self, key, cmd, cwd, answers=None, append=False):
        """Run one step as a process, its output going to logs/<key>.log (added at its end with append)."""
        os.makedirs(cwd, exist_ok=True)
        log_path = os.path.join(self.logdir, key + '.log')
        # The command as written in the log: python and the scripts of this checkout without their full path
        shown = ['python' if c == sys.executable else os.path.relpath(c, REPO_DIR) if c.startswith(REPO_DIR + os.sep) else c
                 for c in cmd]
        with open(log_path, 'a' if append else 'w') as log:
            log.write('$ ' + ' '.join(shown) + '\n')
            if answers is not None:
                log.write('# answers on standard input: {0!r}\n'.format(answers))
            log.flush()
            proc = subprocess.Popen(cmd, cwd=cwd, env=self.env, stdout=log, stderr=subprocess.STDOUT,
                                    stdin=subprocess.PIPE if answers is not None else subprocess.DEVNULL, text=True)
            proc.communicate(answers)
        if proc.returncode == -9:  # SIGKILL: on the server, the out-of-memory killer (memory limit of the container)
            raise StepFailed('{0} was stopped, most likely because it ran out of memory. Please use a more specific '
                             'keyword, or leave some structures out. See logs/{1}.log:\n{2}'.format(
                                 self.status[key]['title'], key, tail(log_path)))
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
        update = self.update_seqfile or seqfile_is_stale(self.seqfile, self.seqres_max_age_days)
        if self.max_entries is not None:
            # Before anything is downloaded: stop if the keyword matches too many PDB entries
            cmd = [sys.executable, os.path.abspath(__file__), '_count', self.keyword, self.rel(self.seqfile),
                   '--max-entries', str(self.max_entries)] + (['--update'] if update else [])
            self.run('download', cmd, cwd=self.projdir)
            update = False
        cmd = [sys.executable, os.path.join(SCRIPTS_DIR, 'CIFcurator_Download.py'), self.keyword,
               '--projdir', '.', '--seqfile', self.rel(self.seqfile)]
        if self.csm:
            cmd += ['--csm', 'True']
            if self.max_models is not None:
                cmd += ['--max-models', str(self.max_models)]
        if update:
            cmd.append('--update')
        self.run('download', cmd, cwd=self.projdir, append=self.max_entries is not None)
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
        cmd = [sys.executable, os.path.abspath(__file__), '_molid', 'clean_bank', 'standard_MolID_bank', '.',
               self.rel(self.seqfile)]
        if self.input_file:
            cmd += ['--input-file', self.rel(self.input_file)]
        if self.excluded_files:
            path = os.path.join(self.logdir, 'excluded_files.txt')
            with open(path, 'w') as f:
                f.write('\n'.join(self.excluded_files) + '\n')
            cmd += ['--exclude-file', self.rel(path)]
        self.run('molid', cmd, cwd=self.projdir)

    def chainid(self):
        os.makedirs(os.path.join(self.projdir, 'standard_ChainID_bank'), exist_ok=True)
        if self.reference_json:
            answers = '2\n5\n{0}\n4\n4\n'.format(self.rel(self.reference_json))
        else:
            answers = '2\n4\n4\n4\n'
        cmd = [sys.executable, os.path.join(SCRIPTS_DIR, 'CIFcurator_ChainStandardization.py'),
               'standard_MolID_bank', 'standard_ChainID_bank']
        self.run('chainid', cmd, cwd=self.projdir, answers=answers)

    def residueid(self):
        os.makedirs(os.path.join(self.projdir, 'standard_ResidueID_bank'), exist_ok=True)
        answers = '1\n4\n{0}\n4\n3\nQUIT\n'.format(self.threshold)
        cwd = os.path.join(self.projdir, 'residue_alignments')
        cmd = [sys.executable, os.path.join(SCRIPTS_DIR, 'CIFcurator_ResidueStandardization.py'),
               self.rel(os.path.join(self.projdir, 'standard_ChainID_bank'), cwd),
               self.rel(os.path.join(self.projdir, 'standard_ResidueID_bank'), cwd)]
        self.run('residueid', cmd, cwd=cwd, answers=answers)
        if not glob.glob(os.path.join(self.projdir, 'standard_ResidueID_bank', '*.cif')):
            raise StepFailed('Step 5 wrote no CIF file; see logs/residueid.log:\n'
                             + tail(os.path.join(self.logdir, 'residueid.log')))

    def reuse(self, key):
        """A step before start_step: its output is already in projdir."""
        if key == 'download':
            if not glob.glob(os.path.join(self.projdir, 'clean_bank', '*.cif')):
                raise StepFailed('Cannot rerun from a later step: there are no downloaded structures (clean_bank) to reuse.')
        elif key == 'exclude':
            path = os.path.join(self.logdir, 'excluded_files.txt')
            if not os.path.isfile(os.path.join(self.projdir, 'excluded_structures.txt')):
                return self.set_state(key, 'skipped')
            if os.path.isfile(path):
                with open(path) as f:
                    self.excluded_files = read_name_list(f.read())
        elif not all(os.path.exists(os.path.join(self.projdir, out)) for out in STEP_OUTPUTS[key] if out.startswith('standard_')):
            raise StepFailed('Cannot rerun from {0}: the output of "{1}" is missing.'.format(
                self.status[self.start_step]['title'], self.status[key]['title']))
        self.set_state(key, 'reused', summary='output of the earlier job')

    def run_all(self, zip_path=None):
        os.makedirs(self.logdir, exist_ok=True)
        with open(os.path.join(self.logdir, 'parameters.json'), 'w') as f:
            json.dump({'keyword': self.keyword, 'threshold': self.threshold, 'start_step': self.start_step,
                       'csm': self.csm, 'max_models': self.max_models, 'exclude': self.exclude,
                       'max_entries': self.max_entries,
                       # relative to projdir (see the top of this file)
                       'input_file': self.input_file and self.rel(self.input_file),
                       'reference_json': self.reference_json and self.rel(self.reference_json),
                       'seqfile': self.rel(self.seqfile),
                       'cifcurator_commit': git_commit()}, f, indent=2)
        self.write_status()
        functions = {'download': self.download, 'exclude': self.exclude_structures, 'molid': self.molid,
                     'chainid': self.chainid, 'residueid': self.residueid}
        keys = [k for k, _ in STEPS if k != 'zip']
        first = keys.index(self.start_step) if self.start_step else 0
        try:
            for i, key in enumerate(keys):
                if i < first:
                    self.reuse(key)
                elif key == 'exclude' and not self.exclude:
                    self.set_state('exclude', 'skipped')
                else:
                    self.step(key, functions[key])
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


# -- Limit on the number of PDB entries ---------------------------------------------------
def count_entries(keyword, seqfile, update=False):
    """
    Number of PDB entries that Step 1 (CIFcurator_Download.py) would download for keyword, found the
    same way as pdbutils.retrieve_sequence_from_PDB(mode='metadata') then
    retrieve_sequences_matching_from_PDB() and get_idset_from_metadata(): the sequences that follow a
    line matching the keyword, then the entries of every chain whose sequence ends with one of them.
    Unlike those functions, it keeps only sets (no numpy arrays of strings, whose memory is the
    longest string times the number of hits: more than 16 GB for "ribosom"). With update, the
    sequence file is downloaded first.
    """
    from CIFcurator import pdbutils
    if update:
        seqfile = pdbutils.retrieve_seqfile(seqfile=seqfile)
    pattern = re.compile(re.escape(keyword), flags=re.I)
    refs_by_length = {}
    with open(seqfile) as f:
        nextline = False
        for line in f:
            if nextline:
                seq = line.strip().upper()
                if seq:
                    refs_by_length.setdefault(len(seq), set()).add(seq)
                nextline = False
            elif pattern.search(line):
                nextline = True
    entries = set()
    with open(seqfile) as f:
        prevline = '#'
        for line in f:
            if not line.startswith('>'):
                seq = line.rstrip('\n').upper()
                for length, refs in refs_by_length.items():
                    if len(seq) >= length and seq[len(seq)-length:] in refs:
                        entries.add(pdbutils.extended_pdbid(prevline.lstrip('>').split()[0].rsplit('_', 1)[0]))
                        break
            prevline = line
    return len(entries)


def check_entries(keyword, seqfile, max_entries, update=False):
    """Exit with a message for the user if keyword matches more than max_entries PDB entries."""
    if update:
        print('Downloading the PDB sequences to {0} ...'.format(seqfile), flush=True)
    n = count_entries(keyword, seqfile, update)
    print('The keyword {0!r} matches {1} PDB entries (this server takes at most {2}).'.format(keyword, n, max_entries))
    if n > max_entries:
        sys.exit('The keyword {0!r} matches {1:,} PDB entries, more than the {2:,} this server can curate in one job. '
                 'Please use a more specific keyword, or run CIFcurator pipeline locally.'.format(keyword, n, max_entries))


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
def link_tree(source, target):
    """Hard links of the file or directory source at target (a copy where hard links are not possible)."""
    if os.path.isfile(source):
        os.makedirs(os.path.dirname(target), exist_ok=True)
        try:
            os.link(source, target)
        except FileExistsError:
            raise
        except OSError:  # e.g. another file system
            shutil.copy2(source, target)
        return
    for root, dirs, files in os.walk(source):
        here = os.path.join(target, os.path.relpath(root, source))
        os.makedirs(here, exist_ok=True)
        for name in files:
            path = os.path.join(root, name)
            if os.path.islink(path):
                os.symlink(os.readlink(path), os.path.join(here, name))
            else:
                link_tree(path, os.path.join(here, name))


def reuse_outputs(source_projdir, target_projdir, start_step):
    """
    Hard links, in target_projdir, of what the steps before start_step wrote in source_projdir
    (STEP_OUTPUTS), so that a rerun from start_step can use them.
    """
    keys = [k for k, _ in STEPS]
    for key in keys[:keys.index(start_step)]:
        for out in STEP_OUTPUTS[key]:
            path = os.path.join(source_projdir, out)
            if os.path.exists(path):
                link_tree(path, os.path.join(target_projdir, out))
    os.makedirs(os.path.join(target_projdir, 'logs'), exist_ok=True)


def step_states(projdir):
    """{step key: state} from projdir/logs/steps.json ({} if there is none)."""
    try:
        with open(os.path.join(projdir, 'logs', 'steps.json')) as f:
            return {s['key']: s['state'] for s in json.load(f)['steps']}
    except (OSError, ValueError, KeyError, TypeError):
        return {}


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
    # A Docker image has no .git folder: deploy/Dockerfile sets CIFCURATOR_COMMIT when it is built
    if os.environ.get('CIFCURATOR_COMMIT'):
        return os.environ['CIFCURATOR_COMMIT']
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

    if argv and argv[0] == '_count':
        p = argparse.ArgumentParser(prog='pipeline.py _count')
        p.add_argument('keyword'); p.add_argument('seqfile')
        p.add_argument('--max-entries', type=int, required=True)
        p.add_argument('--update', action='store_true')
        a = p.parse_args(argv[1:])
        check_entries(a.keyword, a.seqfile, a.max_entries, a.update)
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
    parser.add_argument('--start-step', choices=RERUN_STEPS, default=None,
                        help='Rerun from this step, using the output of the earlier steps already in projdir (or in '
                             '--reuse-from). exclude = leaving structures out, molid = Step 3, chainid = Step 4, '
                             'residueid = Step 5. Default: run every step')
    parser.add_argument('--reuse-from', default=None, metavar='PROJDIR',
                        help='With --start-step: take the output of the earlier steps from this project directory of an '
                             'earlier run, by hard links (it is not changed)')
    parser.add_argument('--max-entries', type=int, default=None, metavar='N',
                        help='Step 1: stop before downloading anything if the keyword matches more than N PDB entries '
                             '(default: no limit)')
    parser.add_argument('--zip', default=None, help='Also write every file of projdir to this zip file')
    a = parser.parse_args(argv)
    if not 1 <= a.threshold <= 100:
        parser.error('--threshold must be between 1 and 100')
    if a.max_models is not None and a.max_models < 1:
        parser.error('--max-models must be at least 1')
    if a.max_entries is not None and a.max_entries < 1:
        parser.error('--max-entries must be at least 1')
    if a.reuse_from and not a.start_step:
        parser.error('--reuse-from needs --start-step')
    exclude = []
    if a.exclude_file:
        with open(a.exclude_file) as f:
            exclude = read_name_list(f.read())
    if a.reuse_from:
        if step_states(a.reuse_from).get('download') not in ('done', 'reused'):
            parser.error('the download step of {0} did not finish, so there is nothing to reuse'.format(a.reuse_from))
        reuse_outputs(os.path.abspath(a.reuse_from), os.path.abspath(a.projdir), a.start_step)
    seqfile = a.seqfile or os.path.join(SERVER_DIR, 'server_data', 'pdb_seqres.fasta')
    os.makedirs(os.path.dirname(os.path.abspath(seqfile)), exist_ok=True)
    pipeline = Pipeline(a.keyword, a.projdir, seqfile, input_file=a.input_file, reference_json=a.reference_json,
                        threshold=a.threshold, csm=a.csm, max_models=a.max_models, exclude=exclude, update_seqfile=a.update, seqres_max_age_days=a.seqres_max_age_days,
                        start_step=a.start_step, max_entries=a.max_entries)
    try:
        pipeline.run_all(zip_path=a.zip)
    except Exception as err:
        print(err, file=sys.stderr)
        return 1
    print('Done. Results in {0}'.format(pipeline.projdir) + (' and {0}'.format(a.zip) if a.zip else ''))
    return 0


if __name__ == '__main__':
    sys.exit(main())
