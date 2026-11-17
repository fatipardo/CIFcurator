#!/Applications/anaconda3/envs/CIFcurator/bin/python
# coding: utf-8
#
# ! ! ! master_molID_class_list is very important
# This is the list that contains every file's MolID class
# !! molIDConversion_list is important . . . it contains the objects
# MolIDConversion and is what is going to be updated by the user and evaulated
# to determine when the next step in the program is unlocked
# Create list of MolIDConversion objects using unique_molID_occur_map

r"""Step 3 of CIFcurator: assign the chain IDs of the entities (MolIDs) found in the CIF files.

Reads every .cif file in source_dir, lists the entities (MolIDs) found in them, and asks
you to give each entity a new chain ID for every copy of it. The rewritten CIF files are
saved in target_dir. The script is interactive: it shows a menu and asks for options.

Conversion Build Menu:
  1) Show full conversion
         Show every entity and the chain IDs assigned to it so far.
  2) Show only unassigned conversions
         Show only the entities that still miss chain IDs.
  3) Enter input file
         Read the chain IDs from a file with one entity per line, in the format
         ENTITY NAME:chainID,chainID,...   for example   WATER:F,G
  4) Generate input file
         From the clustering information, an input file is generated for the user
  5) Continue to next step of curation
         Shown once every entity has all its chain IDs.
         and writes Write OldChainID_NewChainID_numbered_MolID_Files_Map.csv in target_dir
         (original_chain_names:new_chain_names:entity_name:number_of_files:file_name).
  QUIT   Exit without writing CIF files.

After option 5, if no two entities share a chain ID, option 6) Finalize Curation writes the
CIF files to target_dir. Otherwise the Concatenations Menu appears: entities with the same
chain ID will be concatenated into one chain, and each concatenation must be accepted first:
  1) Show all conversions
  2) Show only unaccepted concatenations
  3) Search and modify destination chainIDs of proposed concatenations
  4) Accept proposed concatenation one by one
  5) Accept ALL (make sure you agree with all concatenations)
  6) Finalize Curation

The answers can also be given on standard input, for example:
  printf '3\nTIM_input_file.txt\n7\n6\n' | CIFcurator_MolID.py clean_bank standard_MolID_bank
"""

from __future__ import print_function
import argparse, glob
from CIFcurator import cifcuratormolidcifutils as molidutils
from CIFcurator import cifcluster


########################
# READ INPUT ARGUMENTS #
########################
parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0],
                                 epilog=__doc__.split('\n\n', 1)[1],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('source_dir', help='Directory with the input CIF files (for example, the output of Step 1, CIFcurator_Download.py)')
parser.add_argument('target_dir', help='Directory where the new CIF files are saved')
parser.add_argument('projdir', nargs='?', default='.',
                    help='Project directory, where clustering_storage/ and input_file.txt are written (default: current directory)')
parser.add_argument('seqfile', nargs='?', default=None,
                    help='File with the sequences of all the structures in the PDB (default: {projdir}/seqres.txt)')
args = parser.parse_args()
source_dir = args.source_dir
target_dir = args.target_dir
projdir = args.projdir
seqfile = args.seqfile if args.seqfile is not None else projdir+'/seqres.txt'


#########################################
# READ PDB FILES AND DEFINE MolID LISTS #
#########################################

filelist=glob.glob(source_dir+'/*.cif')
cifcluster.run_clustering_pipeline(filelist, projdir, seqfile)
cluster_pdbchids, outlier_case_entries = cifcluster.finding_problematic_cluster_results(f"{projdir}/clustering_storage/")
master_molID_class_list = molidutils.pdb_to_masterlist(filelist)
unique_molID_occur_map  = molidutils.CreateMasterUniqueMolIDMap(master_molID_class_list)
molIDConversion_list    = molidutils.uniquelist_to_conversionlist(unique_molID_occur_map)
#FAPA MARCH 2024
MolID_to_files_map = molidutils.CreateMasterUniqueMolIDMapWithFileName(master_molID_class_list)
MolID_occur_dict_of_lists = molidutils.CreateMasterUniqueMolIDOccursLIST(master_molID_class_list)
MolID_ChainID_dict_of_lists = molidutils.CreateMasterUniqueMolIDinitialChainIDsLIST(master_molID_class_list)

#####################################
# INTERACTIVE MOLID CONVERSION MENU #
#####################################
input_menu = ""
input_menu_complete = ""
final_menu=""

while(input_menu != "QUIT"):
    if (input_menu_complete == "1"):
        print("""Congratulations! You have successfully constructed your
                 conversion templates. You can proceed to the next section
                 by selection option 7 or, continue to edit your conversion
                 template through this menu
              """)
    print("""CIFcurator MolID Conversion Build Menu
             Select one of the following options to proceed:
             1) Show full conversion
             2) Show only unassigned conversions
             3) Enter input file
             4) Generate input file
             OR type 'QUIT' to exit
          """)
    if (input_menu_complete == "1"):
        print("    5) Continue to next step of curation")
    input_menu = input('Option Number: ')
    if (input_menu == "1"):
        molidutils.show_full_conversion(molIDConversion_list)
    elif (input_menu == "2"):
        molidutils.show_unassigned_conversion(molIDConversion_list)
    elif (input_menu == "3"):
        molIDConversion_list = molidutils.add_user_conversion(molIDConversion_list)
    elif (input_menu == "4"):
        final_chain_set = cifcluster.generate_chain_assignments(molIDConversion_list, master_molID_class_list,
                                                                cluster_pdbchids, outlier_case_entries, projdir)
    elif (input_menu == "5"):
        if (input_menu_complete == "1"):
            master_molID_class_list = molidutils.update_masterlist(master_molID_class_list, molIDConversion_list)
            count_problems = molidutils.problem_counter(master_molID_class_list)
            if (count_problems != 0):
                master_molID_class_list = cifcluster.resolve_polymeric_chain_concatenations(master_molID_class_list, f"{projdir}/updated_input_file.txt")

                # Added additional check
                for molID_class in master_molID_class_list:
                    for ent_id, ent_type in molID_class.entity_type_map.items():
                        if ent_type == 'polymer':
                            ent_name = molID_class.entity_name_map[ent_id]
                            for chID in molID_class.molID_chID.get(ent_name, []):
                                if not molID_class.complete_order.get(chID, True):
                                    raise RuntimeError(
                                        f"Unresolved polymer chain collision in {molID_class.file_name}: "
                                        f"chain {chID} ({ent_name}) -> {molID_class.chID_newchID_map[chID]}"
                                    )

                allnewchains = molidutils.return_unassigned_conversion(master_molID_class_list, step='concatenation')
                for newchain in allnewchains:
                    master_molID_class_list = molidutils.list_accept_concatenations_auto(
                        master_molID_class_list, newchain, new_order=None)[0]

            print("Saving 'Track changes' to output directory")
            molidutils.show_full_conversion_and_file_list_by_number_chains(
                master_molID_class_list, MolID_to_files_map, target_dir
            )

            input_menu = "QUIT"
            final_menu = "START"
    input_menu_complete = molidutils.check_complete(molIDConversion_list)

if (final_menu == "START"):

    count_problems = molidutils.problem_counter(master_molID_class_list)
    if (count_problems == 0):
        final_menu_complete = "1"

    if (final_menu_complete == "1"):
        print("    6) Finalize Curation")

    final_menu = input('Option Number: ')

    if (final_menu == "6"):
        print("Finalizing Curation ...")
        molidutils.masterlist_to_pdb(filelist, master_molID_class_list, target_dir=target_dir)
        final_menu = "QUIT"
    else:
        print("Sorry, something went wrong, try again")