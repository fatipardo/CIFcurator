#!/usr/bin/env python
# coding: utf-8
#

r"""Step 1 of CIFcurator: download the structures that match a keyword from the PDB, and clean them.

It does two things, without asking for input (they used to be Steps 1 and 2, run from notebooks):

Download structural ensemble from RCSB PDB
  - Downloads the sequences of all the structures in the PDB (pdb_seqres.fasta.gz, from the
    wwPDB beta archive) to {projdir}/seqres.txt. If that file already exists it is reused;
    use --update to download it again.
  - Finds the reference sequences whose metadata matches the keyword, and then all the chains
    in the PDB whose sequence matches one of them.
  - Downloads the biological assemblies of those structures to {projdir}/raw_bank, one file
    per assembly, named PDBID+NN.cif (e.g. pdb_00001aw1+01.cif).
  - Optional (--csm True): downloads the computed structure models (CSMs) that match the
    keyword, from AlphaFold DB, to {projdir}/raw_bank, named AFID+01.cif
    (e.g. AF-P60174-F1+01.cif). A CSM contains ONE PROTEIN CHAIN ONLY: no nucleic acids,
    ligands, waters or other copies of the protein.

Clean structures: one CIF file per biological assembly
  - Rewrites every CIF file of {projdir}/raw_bank to {projdir}/clean_bank, with only a
    limited set of categories, written in a standard way.

The output of this script, {projdir}/clean_bank, is the input of Step 3 (CIFcurator_MolID.py).

Examples (use quotes when the keyword has more than one word):
  CIFcurator_Download.py "triosephosphate isomerase" --projdir ./TIM
  CIFcurator_Download.py "triosephosphate isomerase" --projdir ./TIM --csm True --max-models 10
  CIFcurator_Download.py "triose[- ]?phosphate isomerase" --projdir ./TIM --regex
  CIFcurator_Download.py "h(a)?emoglobin" --projdir ./HB --regex --csm True --csm-keyword "hemoglobin"
"""

import argparse, os, re
from CIFcurator import cifcurator_io, pdbutils, cleanutils


def str_to_bool(value):
    """Reads True/False (any case) from the command line."""
    if value.lower() == 'true':
        return True
    if value.lower() == 'false':
        return False
    raise argparse.ArgumentTypeError('expected True or False, got {0!r}'.format(value))


########################
# READ INPUT ARGUMENTS #
########################
parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0],
                                 epilog=__doc__.split('\n\n', 1)[1],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('keyword', help='Keyword to search for in the PDB, e.g. "triosephosphate isomerase" (in quotes)')
parser.add_argument('--projdir', default='.',
                    help='Project directory, where raw_bank and clean_bank are created (default: current directory)')
parser.add_argument('--csm', type=str_to_bool, default=False, metavar='{True,False}',
                    help='Also download computed structure models (CSMs) from AlphaFold DB (default: False)')
parser.add_argument('--max-models', type=int, default=None, metavar='N',
                    help='With --csm True, download at most N CSMs, best matches first (default: all of them)')
parser.add_argument('--csm-keyword', default=None, metavar='TEXT',
                    help='With --csm True, plain text for the CSM search, instead of the keyword '
                         '(required with --regex)')
parser.add_argument('--seqfile', default=None,
                    help='File with the sequences of all the structures in the PDB (default: {projdir}/seqres.txt)')
parser.add_argument('--update', action='store_true',
                    help='Download the sequence file again, even if it already exists')
parser.add_argument('--regex', action='store_true',
                    help='Advanced: the keyword is a Python regular expression, not plain text '
                         '(with --csm True, also give --csm-keyword)')
args = parser.parse_args()
if args.regex:
    if args.csm and args.csm_keyword is None:
        parser.error('--regex with --csm True needs --csm-keyword: the CSM search of RCSB PDB does not '
                     'understand regular expressions, so give it plain text, e.g. --csm-keyword "hemoglobin"')
    try:
        re.compile(args.keyword)
    except re.error as err:
        parser.error('the keyword {0!r} is not a valid regular expression: {1}'.format(args.keyword, err))
keyword = args.keyword
projdir = args.projdir
seqfile = args.seqfile if args.seqfile is not None else projdir+'/seqres.txt'


################################
# DOWNLOAD STRUCTURES FROM PDB #
################################

cifcurator_io.check_project(projdir=projdir)

update = args.update or not os.path.exists(seqfile)
if update:
    print('Downloading the PDB sequences to {0} ...'.format(seqfile))
else:
    print('Using the PDB sequences in {0} (use --update to download them again)'.format(seqfile))

ref_sequences, ref_metadata = pdbutils.retrieve_sequence_from_PDB(keyword, mode='metadata', update=update, seqfile=seqfile,
                                                                  regex=args.regex)
print('{0} sequences were identified as potential hits!'.format(len(ref_sequences)))
if len(ref_sequences) == 0 and not args.csm:
    parser.exit(1, 'No sequence in the PDB matches the keyword {0!r}. Nothing to download.\n'.format(keyword))

sequences, metadata = pdbutils.retrieve_sequences_matching_from_PDB(ref_sequences, seqfile)
print('{0} sequences were retrieved!'.format(len(sequences)))
for iseq in range(len(sequences)):
    print('{0} {1}'.format(metadata[iseq].strip(), sequences[iseq].strip()))

cifcurator_io.check_project(projdir=projdir, level='raw_bank')
pdbids = pdbutils.get_idset_from_metadata(metadata)
print('Downloading the biological assemblies of {0} PDB entries to {1}/raw_bank ...'.format(len(pdbids), projdir))
pdbutils.download_pdb_from_metadata(metadata, projdir=projdir)

if args.csm:
    csm_keyword = args.csm_keyword if args.csm_keyword is not None else keyword
    pdbutils.download_csm_from_keyword(csm_keyword, projdir=projdir, max_models=args.max_models)
else:
    for option, value in (('--max-models', args.max_models), ('--csm-keyword', args.csm_keyword)):
        if value is not None:
            print('{0} is ignored, because --csm is False'.format(option))


###############################################
# CLEAN STRUCTURES: ONE CIF FILE PER ASSEMBLY #
###############################################

cifcurator_io.check_project(projdir=projdir, level='clean_bank')
print('Cleaning the CIF files of {0}/raw_bank to {0}/clean_bank ...'.format(projdir))
cleanutils.process(projdir=projdir, step='clean', source='raw_bank', target='clean_bank')
print('Done. The clean CIF files are in {0}/clean_bank'.format(projdir))
