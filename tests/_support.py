"""
Shared helpers for the CIFcurator test suite.

Importing this module makes the repository's src/ importable as the CIFcurator
package (so the tests always exercise this checkout, not an installed copy).

Environment variables
---------------------
CIFCURATOR_SKIP_NETWORK=1         skip the pipeline tests (they download from the PDB)
CIFCURATOR_TEST_CACHE=<dir>       where pdb_seqres.txt is cached (default: tests/.cache)
CIFCURATOR_SEQRES_MAX_AGE_DAYS=N  re-download pdb_seqres.txt when older than N days (default: 30)
CIFCURATOR_KEEP_TEST_OUTPUT=1     keep the pipeline work directory after the run
CIFCURATOR_SCRIPT_TIMEOUT=<sec>   timeout for each interactive script (default: 1800)
"""
import atexit
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
    # cleanutils.simplify_cif never closes the files it writes (CPython closes them on reassignment).
    warnings.filterwarnings('ignore', category=ResourceWarning)


def seqres_file():
    """
    Return the path to a cached pdb_seqres.txt, downloading it the first time
    or when the cached copy is older than SEQRES_MAX_AGE_DAYS.

    The download goes to a temporary name and is renamed only once complete, so an
    interrupted download never leaves a truncated file that later runs would trust.
    retrieve_seqfile raises ContentTooShortError when fewer bytes arrive than the
    server's Content-Length, so a truncated response is never cached either.
    """
    from CIFcurator import pdbutils
    os.makedirs(CACHE_DIR, exist_ok=True)
    path = os.path.join(CACHE_DIR, 'pdb_seqres.txt')
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
            raise RuntimeError('Downloaded pdb_seqres.txt does not look like FASTA: {0!r}'.format(first[:80]))
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
            os.killpg(proc.pid, signal.SIGKILL)
            proc.wait()
            problem = 'timed out after {0:.0f} s'.format(timeout)
    if problem:
        raise ScriptFailed('{0} {1}. Last lines of {2}:\n{3}'.format(
            script, problem, log_path, tail(log_path)))


def tail(path, n=40):
    with open(path, errors='replace') as f:
        return ''.join(f.readlines()[-n:])


def cif_columns(path, *keys):
    """Return the requested mmCIF items of `path` as lists (via Biopython)."""
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict
    d = MMCIF2Dict(path)
    return [d[k] if isinstance(d[k], list) else [d[k]] for k in keys]
