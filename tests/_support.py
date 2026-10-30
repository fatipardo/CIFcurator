"""
Shared helpers for the CIFcurator test suite.

Importing this module makes the repository's src/ importable as the CIFcurator
package (so the tests always exercise this checkout, not an installed copy).

Environment variables
---------------------
CIFCURATOR_SKIP_NETWORK=1         skip the pipeline tests (they download from the PDB)
CIFCURATOR_TEST_CACHE=<dir>       where pdb_seqres.fasta is cached (default: tests/.cache)
CIFCURATOR_SEQRES_MAX_AGE_DAYS=N  re-download pdb_seqres.fasta when older than N days (default: 30)
CIFCURATOR_KEEP_TEST_OUTPUT=1     keep the pipeline work directory after the run
CIFCURATOR_SCRIPT_TIMEOUT=<sec>   timeout for each interactive script (default: 1800)
"""
import atexit
import contextlib
import os
import shutil
import signal
import subprocess
import sys
import tempfile
import time
import warnings

TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_DIR = os.path.dirname(TESTS_DIR)
SRC_DIR = os.path.join(REPO_DIR, 'src')
SCRIPTS_DIR = os.path.join(REPO_DIR, 'scripts')
DATA_DIR = os.path.join(TESTS_DIR, 'data')

CACHE_DIR = os.environ.get('CIFCURATOR_TEST_CACHE', os.path.join(TESTS_DIR, '.cache'))
SEQRES_MAX_AGE_DAYS = float(os.environ.get('CIFCURATOR_SEQRES_MAX_AGE_DAYS', '30'))
SKIP_NETWORK = os.environ.get('CIFCURATOR_SKIP_NETWORK', '') not in ('', '0')
KEEP_OUTPUT = os.environ.get('CIFCURATOR_KEEP_TEST_OUTPUT', '') not in ('', '0')
SCRIPT_TIMEOUT = float(os.environ.get('CIFCURATOR_SCRIPT_TIMEOUT', '1800'))

# The library has many non-raw regex strings; keep test output readable.
warnings.filterwarnings('ignore', category=SyntaxWarning)

# The package lives in src/ but is imported as CIFcurator (see setup.py package_dir).
# Expose it under that name through a symlink in a private directory, used both
# in-process (sys.path) and by the scripts run as subprocesses (PYTHONPATH).
_SHIM_DIR = tempfile.mkdtemp(prefix='cifcurator_tests_pylib_')
os.symlink(SRC_DIR, os.path.join(_SHIM_DIR, 'CIFcurator'))
sys.path.insert(0, _SHIM_DIR)
atexit.register(shutil.rmtree, _SHIM_DIR, ignore_errors=True)


def quiet_known_warnings():
    """
    Hide warnings about known, harmless issues in the library, so real failures stand out.
    Call from setUpModule: unittest resets warning filters when the run starts.
    """
    warnings.filterwarnings('ignore', category=SyntaxWarning)
    # Library code that leaves files for CPython to close.
    warnings.filterwarnings('ignore', category=ResourceWarning)


def seqres_file():
    """
    Return the path to a cached pdb_seqres.fasta (wwPDB beta archive, extended PDB IDs), downloading it the first time
    or when the cached copy is older than SEQRES_MAX_AGE_DAYS.

    The download goes to a temporary name and is renamed only once complete, so an
    interrupted download never leaves a truncated file that later runs would trust.
    retrieve_seqfile raises ContentTooShortError when fewer bytes arrive than the
    server's Content-Length, so a truncated response is never cached either.
    """
    from CIFcurator import pdbutils
    os.makedirs(CACHE_DIR, exist_ok=True)
    path = os.path.join(CACHE_DIR, 'pdb_seqres.fasta')  # OCT2026: not .txt, so a legacy cache is not reused
    if os.path.isfile(path):
        age_days = (time.time() - os.path.getmtime(path)) / 86400
        if age_days < SEQRES_MAX_AGE_DAYS:
            return path
    fd, partial = tempfile.mkstemp(prefix='pdb_seqres.', suffix='.part', dir=CACHE_DIR)
    os.close(fd)
    try:
        pdbutils.retrieve_seqfile(seqfile=partial)
        with open(partial) as f:
            first = f.readline()
        if not first.startswith('>'):
            raise RuntimeError('Downloaded pdb_seqres.fasta does not look like FASTA: {0!r}'.format(first[:80]))
        os.replace(partial, path)
    finally:
        if os.path.exists(partial):
            os.remove(partial)
    return path


class ScriptFailed(AssertionError):
    pass


def run_script(script, args, answers, cwd, log_path, timeout=None):
    """
    Run one of the interactive scripts in scripts/, feeding `answers` on stdin.

    The script runs with `cwd` as its working directory: the alignment code writes
    its MUSCLE input/output files (Seq.fa, <chain>.fa, <chain>.fasta) in the current
    directory and reuses an existing <chain>.fasta, so every step needs its own
    empty directory.

    Output (stdout+stderr) is written to `log_path`. Raises ScriptFailed with the end
    of the log if the script exits non-zero or exceeds the timeout. (If MUSCLE fails,
    AlignSequences_v4 waits forever for its output, so the timeout matters.)
    """
    timeout = SCRIPT_TIMEOUT if timeout is None else timeout
    os.makedirs(cwd, exist_ok=True)
    env = os.environ.copy()
    env['PYTHONPATH'] = os.pathsep.join(p for p in (_SHIM_DIR, env.get('PYTHONPATH')) if p)
    env['PYTHONWARNINGS'] = 'ignore::SyntaxWarning'
    cmd = [sys.executable, os.path.join(SCRIPTS_DIR, script)] + list(args)
    with open(log_path, 'w') as log:
        # New session, so a timeout also kills MUSCLE started by the script.
        proc = subprocess.Popen(cmd, cwd=cwd, env=env, stdin=subprocess.PIPE, stdout=log,
                                stderr=subprocess.STDOUT, text=True, start_new_session=True)
        try:
            proc.communicate(answers, timeout=timeout)
            problem = None if proc.returncode == 0 else 'exited with code {0}'.format(proc.returncode)
        except subprocess.TimeoutExpired:
            _kill_group(proc)
            problem = 'timed out after {0:.0f} s'.format(timeout)
        except BaseException:
            # Ctrl-C reaches only the test runner (the script is in its own session),
            # so stop the script and MUSCLE here too before passing the exception on.
            _kill_group(proc)
            raise
    if problem:
        raise ScriptFailed('{0} {1}. Last lines of {2}:\n{3}'.format(
            script, problem, log_path, tail(log_path)))


def _kill_group(proc):
    """Kill the process group started by run_script (the script and MUSCLE) and reap the script."""
    with contextlib.suppress(ProcessLookupError):
        os.killpg(proc.pid, signal.SIGKILL)
    proc.wait()


def tail(path, n=40):
    with open(path, errors='replace') as f:
        return ''.join(f.readlines()[-n:])


def cif_columns(path, *keys):
    """Return the requested mmCIF items of `path` as lists (via Biopython)."""
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict
    d = MMCIF2Dict(path)
    return [d[k] if isinstance(d[k], list) else [d[k]] for k in keys]


def mmcif_problems(path):
    """
    Return a list of problems that make `path` inconsistent with the mmCIF dictionary, as pointed out
    by the manuscript reviewer (OCT2026): missing categories that the atoms refer to, missing key items,
    and links between categories that do not match (e.g. a chain in _atom_site with no _struct_asym row,
    or _pdbx_poly_seq_scheme numbering that differs from the atoms after Steps 3-5).
    An empty list means no problem was found.
    """
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict
    d = MMCIF2Dict(path)

    def col(key):
        value = d.get(key)
        return None if value is None else (value if isinstance(value, list) else [value])

    def null(value):
        return value in ('?', '.')

    problems = []
    required = ['_atom_type.symbol', '_chem_comp.id', '_struct_asym.id', '_struct_asym.entity_id',
                '_entity.id', '_entity_poly.entity_id', '_entity_poly_seq.num', '_pdbx_poly_seq_scheme.asym_id',
                '_citation.id', '_citation_author.citation_id', '_citation_author.ordinal',
                '_audit_conform.dict_name', '_audit_conform.dict_version']
    # OCT2026: computed structure models (RCSB IDs AF_... and MA_...) have no experimental method
    if not (col('_entry.id') or [''])[0].upper().startswith(('AF_', 'MA_')):
        required.append('_exptl.entry_id')
    for key in required:
        if col(key) is None:
            problems.append('missing ' + key)
    # OCT2026: the dictionary version the file declares is the one CIFcurator follows
    from CIFcurator import cleanutils
    if col('_audit_conform.dict_version') not in (None, [cleanutils.MMCIF_DICT_VERSION]):
        problems.append('_audit_conform.dict_version is not ' + cleanutils.MMCIF_DICT_VERSION)
    # OCT2026: not a dictionary item
    if col('_exptl.resolution') is not None:
        problems.append('_exptl.resolution is not a dictionary item')
    if problems:
        return problems

    atoms = list(zip(col('_atom_site.type_symbol'), col('_atom_site.label_comp_id'), col('_atom_site.label_asym_id'),
                     col('_atom_site.label_entity_id'), col('_atom_site.label_seq_id'),
                     col('_atom_site.pdbx_PDB_ins_code'), col('_atom_site.auth_seq_id'),
                     col('_atom_site.auth_asym_id')))
    asym_entity = {}
    for asym, entity in zip(col('_struct_asym.id'), col('_struct_asym.entity_id')):
        if asym in asym_entity:
            problems.append('duplicate _struct_asym.id ' + asym)
        asym_entity[asym] = entity
    poly_seq = set(zip(col('_entity_poly_seq.entity_id'), col('_entity_poly_seq.num'), col('_entity_poly_seq.mon_id')))
    scheme = {}
    for row in zip(col('_pdbx_poly_seq_scheme.asym_id'), col('_pdbx_poly_seq_scheme.entity_id'),
                   col('_pdbx_poly_seq_scheme.seq_id'), col('_pdbx_poly_seq_scheme.mon_id'),
                   col('_pdbx_poly_seq_scheme.pdb_seq_num'), col('_pdbx_poly_seq_scheme.pdb_strand_id'),
                   col('_pdbx_poly_seq_scheme.pdb_ins_code')):
        asym, entity, seq, mon, pdb_num, strand, ins = row
        if asym_entity.get(asym) != entity:
            problems.append('_pdbx_poly_seq_scheme chain {0} entity {1} not in _struct_asym'.format(asym, entity))
        if (entity, seq, mon) not in poly_seq:
            problems.append('_pdbx_poly_seq_scheme {0} {1} {2} not in _entity_poly_seq'.format(asym, seq, mon))
        scheme[(asym, seq, mon)] = (pdb_num, strand, '' if null(ins) else ins)

    checks = [('_atom_site.type_symbol', {a[0] for a in atoms}, set(col('_atom_type.symbol'))),
              ('_atom_site.label_comp_id', {a[1] for a in atoms}, set(col('_chem_comp.id'))),
              ('_atom_site.label_asym_id', {a[2] for a in atoms}, set(asym_entity)),
              ('_atom_site.label_entity_id', {a[3] for a in atoms}, set(col('_entity.id'))),
              ('_struct_asym.entity_id', set(asym_entity.values()), set(col('_entity.id'))),
              ('_entity_poly_seq.entity_id', set(col('_entity_poly_seq.entity_id')), set(col('_entity_poly.entity_id'))),
              ('_citation_author.citation_id', set(col('_citation_author.citation_id')), set(col('_citation.id'))),
              ('_exptl.entry_id', set(col('_exptl.entry_id') or []), set(col('_entry.id')))]
    for item, values, parents in checks:
        if values - parents:
            problems.append('{0} without parent: {1}'.format(item, sorted(values - parents)[:5]))

    seen = set()
    for symbol, comp, asym, entity, seq, ins, auth_seq, auth_asym in atoms:
        if (asym, seq, comp) in seen:
            continue
        seen.add((asym, seq, comp))
        if asym in asym_entity and asym_entity[asym] != entity:
            problems.append('chain {0}: atoms in entity {1}, _struct_asym says {2}'.format(asym, entity, asym_entity[asym]))
        if null(seq):
            continue
        if (entity, seq, comp) not in poly_seq:
            problems.append('residue {0} {1} {2} not in _entity_poly_seq'.format(asym, seq, comp))
        expected = (auth_seq, auth_asym, '' if null(ins) else ins)
        if scheme.get((asym, seq, comp)) != expected:
            problems.append('residue {0} {1} {2}: atoms say {3}, _pdbx_poly_seq_scheme says {4}'.format(
                asym, seq, comp, expected, scheme.get((asym, seq, comp))))
    return problems
