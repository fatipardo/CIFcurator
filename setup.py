"""
setup.py: Install CIFcurator
"""

import os
import sys
import re
import subprocess
from os.path import join as pjoin
from glob import glob

from setuptools import setup

# ------------------------------------------------------------------------------
# HEADER
#

VERSION      = "1.0.0"
ISRELEASED   = False
DISABLE_CUDA = True
__author__   = "Levitt Lab, Stanford"
__version__  = VERSION

metadata = {
    'name': 'CIFcurator',
    'version': VERSION,
    'author': __author__,
    'author_email': 'fpardo@stanford.edu',
    'license': 'MIT',
    'url': 'https://github.com/fatipardo/CIFcurator',
    'download_url': 'https://github.com/fatipardo/CIFcurator',
    'platforms': ['Linux', 'OSX'],
    'description': "CIF curation tools",
    'long_description': """CIFcurator offers curation tools for structural ensemble deposited in the Protein Data Bank."""}

# ------------------------------------------------------------------------------
# HELPER FUNCTIONS -- path finding, git, python version, readthedocs
#

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'


def print_warning(string):
    print(bcolors.WARNING + string + bcolors.ENDC)


def find_in_path(name, path):
    "Find a file in a search path"
    #adapted fom http://code.activestate.com/recipes/52224-find-a-file-given-a-search-path/
    for dir in path.split(os.pathsep):
        binpath = pjoin(dir, name)
        if os.path.exists(binpath):
            return os.path.abspath(binpath)
    return None


def git_version():
    """
    Return the git revision as a string.
    Copied from numpy setup.py
    """

    def _minimal_ext_cmd(cmd):
        # construct minimal environment
        env = {}
        for k in ['SYSTEMROOT', 'PATH']:
            v = os.environ.get(k)
            if v is not None:
                env[k] = v
        # LANGUAGE is used on win32
        env['LANGUAGE'] = 'C'
        env['LANG'] = 'C'
        env['LC_ALL'] = 'C'
        out = subprocess.Popen(cmd, stdout = subprocess.PIPE, env=env).communicate()[0]
        return out

    try:
        out = _minimal_ext_cmd(['git', 'rev-parse', 'HEAD'])
        GIT_REVISION = out.strip().decode('ascii')
    except OSError:
        GIT_REVISION = "Unknown"

    return GIT_REVISION

# -----------------------------------------------------------------------------
# INSTALL

metadata['packages']     = ['CIFcurator']
metadata['package_dir']  = {'CIFcurator' :         'src'}
metadata['ext_modules']  = []
metadata['scripts']      = [s for s in glob('scripts/*') if not s.endswith('__.py')]
metadata['python_requires'] = '>=3.10'
#metadata['data_files']   = [('reference', glob('./reference/*'))]
#metadata['cmdclass']     = {'build_ext': custom_build_ext}

# ------------------------------------------------------------------------------
#
# Finally, print a warning at the *end* of the build if something fails
#

def print_warnings():
    print("\n")

if __name__ == '__main__':
    setup(**metadata) # ** will unpack dictionary 'metadata' providing the values as arguments
    print_warnings()
