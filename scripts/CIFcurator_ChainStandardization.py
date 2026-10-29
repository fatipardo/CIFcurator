#!/usr/bin/env python
# coding: utf-8
#
r"""Step 4 of CIFcurator: standardize the chain IDs of the CIF files.

Reads every .cif file in source_dir and builds a set of Standard Sequences, one per chain
ID. The chains of every structure are then aligned (pairwise) to the Standard Sequences and
renamed if needed, so that the same molecule has the same chain ID in every file. The new
CIF files are saved in target_dir. The script is interactive: it shows a menu and asks for
options.

ChainID Standardization Menu:
  1) Select Standard Sequences from a chosen input structure
         Use the chains of one input structure as the Standard Sequences.
  2) Generate Standard Sequences based on all the input structures
         For each chain ID, use the sequence with the highest consensus score among all
         the input structures. Its submenu also has
         5) Load previously generated standard sequences
         to read them from a StandardSequences.json file written in an earlier run.
  3) Inspect/Edit Standard Sequences
         Shown once Standard Sequences exist. Inspect or edit the sequence of a chain ID.
  4) Perform Standardization of Chain IDs
         Shown once Standard Sequences exist. Chain IDs can be added to an ignore list;
         then option 4 of the submenu aligns the chains and writes the CIF files.
  QUIT   Exit without writing CIF files.

The Standard Sequences are saved as StandardSequences.json in target_dir every time they
change. MUSCLE alignment files are written in the current directory.

The answers can also be given on standard input, for example:
  printf '2\n4\n4\n4\n' | CIFcurator_ChainStandardization.py standard_MolID_bank standard_ChainID_bank
"""

from __future__ import print_function
from __future__ import division
import argparse, glob, json
from CIFcurator import cifcuratorchainstandardizationutils as chainstd

########################
# READ INPUT ARGUMENTS #
########################
parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0],
                                 epilog=__doc__.split('\n\n', 1)[1],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('source_dir', help='Directory with the input CIF files (for example, the output of Step 3)')
parser.add_argument('target_dir', help='Directory where the new CIF files and StandardSequences.json are saved')
args = parser.parse_args()
source_dir=args.source_dir
target_dir=args.target_dir


#############################################
# READ PDB FILES AND DEFINE STRUCTURE LISTS #
#############################################

filelist=glob.glob(source_dir+'/*.cif')
Structure_Sequences, structid_list, chid_list = chainstd.pdb_to_structurelists(filelist)
Standard_Sequences = {}


# OCT2026: the standard sequences are saved as JSON in the target directory, so that they
# can be loaded again with option 2) -> 5) without editing the file
def save_standard_sequences(Standard_Sequences):
    if Standard_Sequences:
        std_seq_file = chainstd.save_standard_sequences(Standard_Sequences, target_dir)
        print("Standard sequences saved to " + std_seq_file)


############################################
# INTERACTIVE ChainID STANDARDIZATION MENU #
############################################
input_menu = ""
input_menu_check_1 = ""

while(input_menu != "QUIT"):
    print("CIFcurator ChainID Standardization Menu",
          "    Select one of the following options to proceed:",
          "    1) Select Standard Sequences from a chosen input structure",
          "    2) Generate Standard Sequences based on all the input structures",
          sep="\n")
    if(input_menu_check_1 == "1"):
        print("    3) Inspect/Edit Standard Sequences",
              "    4) Perform Standardization of Chain IDs",
              sep="\n")
    input_menu = input('Option Number: ')
    if (input_menu == "1"):
        Standard_Sequences, input_menu_check_1 = chainstd.select_standard_seq_from_reference(Structure_Sequences,
                                                                                             Standard_Sequences,
                                                                                             structid_list,
                                                                                             input_menu_check_1)
        save_standard_sequences(Standard_Sequences)
    elif (input_menu == "2"):
        Standard_Sequences, input_menu_check_1 = chainstd.create_standard_seq_from_consensus(Structure_Sequences,
                                                                                             Standard_Sequences,
                                                                                             chid_list,
                                                                                             input_menu_check_1)
        print("These are the standard sequences:")
        print(json.dumps(Standard_Sequences, indent=4))
        save_standard_sequences(Standard_Sequences)
    elif (input_menu == "3" and input_menu_check_1 == "1"):
        chainstd.review_standard_seq(Structure_Sequences, Standard_Sequences)
        save_standard_sequences(Standard_Sequences)

    elif (input_menu == "4" and input_menu_check_1=="1"):
        chainstd.align_to_std_seq_and_save_to_disk(Structure_Sequences,
                                                   Standard_Sequences,
                                                   structid_list,
                                                   filelist,
                                                   target_dir=target_dir)
        print("Done!")
        input_menu = "QUIT"
