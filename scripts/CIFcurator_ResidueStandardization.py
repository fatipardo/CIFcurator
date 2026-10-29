#!/Users/fatima/anaconda3/envs/PDBCleanV2/bin/python
# coding: utf-8

r"""Step 5 of CIFcurator: standardize the residue numbers of the CIF files.

Reads every .cif file in source_dir and performs a multiple sequence alignment (MUSCLE) of
each chain ID across all the structures. Residues in the same column of the alignment get
the same residue number in every file. The new CIF files are saved in target_dir. The
script is interactive: it shows a menu and asks for options.

Residue Number Standardization Menu:
  1) Proceed to multiple alignment menu
         Choose the chain IDs to standardize and perform the multiple alignments. You are
         asked for an occupancy threshold (1-100): residue columns present in at most
         that percentage of the structures are numbered with letters (insertion codes).
  2) View conversion template
         Show the old and new residue numbers of every structure.
  3) Perform residue number standardization
         Write the CIF files with the new residue numbers to target_dir.
  4) Save conversion template
         Write OldResID_NewResID_Map.csv in target_dir.
  QUIT   Exit.
Options 2 to 4 are shown after the multiple alignments are done.

MUSCLE alignment files, and residue_warning_structures.csv if some structures need
two-letter insertion codes, are written in the current directory.

The answers can also be given on standard input, for example:
  printf '1\n1\n4\n75\n4\n3\nQUIT\n' | CIFcurator_ResidueStandardization.py standard_ChainID_bank standard_ResidueID_bank
"""

from __future__ import print_function
from __future__ import division
import argparse, glob
from CIFcurator import cifcuratorresiduestandardizationutils as resstd

########################
# READ INPUT ARGUMENTS #
########################
parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0],
                                 epilog=__doc__.split('\n\n', 1)[1],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('source_dir', help='Directory with the input CIF files (for example, the output of Step 4)')
parser.add_argument('target_dir', help='Directory where the new CIF files and OldResID_NewResID_Map.csv are saved')
args = parser.parse_args()
source_dir=args.source_dir
target_dir=args.target_dir


#############################################
# READ PDB FILES AND DEFINE STRUCTURE LISTS #
#############################################
filelist=glob.glob(source_dir+'/*.cif')
Structure_Sequences, ChID_ResiNum_Vector, structid_list, chid_list = resstd.pdb_to_structurelists(filelist)


############################################
# INTERACTIVE RESIDUE STANDARDIZATION MENU #
############################################
input_menu = ""
input_menu_check = ""

while(input_menu != "QUIT"):
    print("CIFcurator Residue Number Standardization Menu",
          "    After checking all structures are loaded, select option 1 to proceed:",
          "    1) Proceed to multiple alignment menu",
          sep="\n")
    if(input_menu_check == "1"):
        print("    2) View conversion template")
        print("    3) Perform residue number standardization")
        print("    4) Save conversion template")
        print("    OR Type QUIT to exit")
    input_menu = input('Option Number: ')
    if (input_menu == "1"):
        Structure_Sequences_Aligned, Structure_ConversionTemplate, chid_list, input_menu_check = resstd.perform_multiple_alignment(Structure_Sequences,
                                                                                                                                   ChID_ResiNum_Vector,
                                                                                                                                   structid_list,
                                                                                                                                   chid_list,
                                                                                                                                   input_menu_check)
    elif (input_menu == "2" and input_menu_check == "1"):
        resstd.show_conversiontemplate(Structure_ConversionTemplate)
    elif (input_menu == "4" and input_menu_check == "1"):
        resstd.write_and_show_conversiontemplate(Structure_ConversionTemplate,target_dir,True)
    elif (input_menu == "3" and input_menu_check == "1"):
        resstd.conversiontemplate_to_pdb_FAPA(filelist, Structure_ConversionTemplate, target_dir=target_dir)
