#!/Users/fatima/anaconda3/envs/PDBCleanV2/bin/python
# coding: utf-8
#
# ! ! ! master_molID_class_list is very important
# This is the list that contains every file's MolID class
# !! molIDConversion_list is important . . . it contains the objects
# MolIDConversion and is what is going to be updated by the user and evaulated
# to determine when the next step in the program is unlocked
# Create list of MolIDConversion objects using unique_molID_occur_map

from __future__ import print_function
import sys, glob
from CIFcurator import cifcuratormolidcifutils as molidutils


########################
# READ INPUT ARGUMENTS #
########################
n_arg = len(sys.argv)
if(n_arg<3):
    print('Usage error: {0} <source directory> <target directory>'.format(sys.argv[0]))
    sys.exit()
source_dir = sys.argv[1]
target_dir = sys.argv[2]


#########################################
# READ PDB FILES AND DEFINE MolID LISTS #
#########################################

filelist=glob.glob(source_dir+'/*.cif')
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
# Goal:
# users complete their molID conversion templates by ensuring that each member of
# molIDConversion_list has status complete = True
input_menu = ""
input_menu_complete = ""
# For use in the next section
concat_menu = ""
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
             4) Search MolID to add chain ID conversion
             5) Go entry by entry to add chain ID conversion
             6) Remove a chain ID conversion
             A) Track changes (original_chain_name:new_chain:entity:file_name)
          """)
    if (input_menu_complete == "1"):
        print("    7) Continue to next step of curation")
    input_menu = input('Option Number: ')
    if (input_menu == "1"):
        molidutils.show_full_conversion(molIDConversion_list)
    elif (input_menu == "2"):
        molidutils.show_unassigned_conversion(molIDConversion_list)
    elif (input_menu == "3"):
        molIDConversion_list = molidutils.add_user_conversion(molIDConversion_list)
    elif (input_menu == "4"):
        molIDConversion_list = molidutils.edit_conversion_interface(molIDConversion_list, action='add')#FAPA
    elif (input_menu == "5"):
        molIDConversion_list = molidutils.edit_conversion_manual(molIDConversion_list)
    elif (input_menu == "6"):
        molIDConversion_list = molidutils.edit_conversion_interface(molIDConversion_list, action='remove')
    elif (input_menu == "A"):
        molidutils.show_full_conversion_and_file_list_by_number_chains(MolID_ChainID_dict_of_lists,
                                                                       molIDConversion_list, MolID_to_files_map,
                                                                       MolID_occur_dict_of_lists, target_dir)
    elif (input_menu == "7"):
        if (input_menu_complete == "1"):
            master_molID_class_list = molidutils.update_masterlist(master_molID_class_list, molIDConversion_list)
            count_problems = molidutils.problem_counter(master_molID_class_list)
            if (count_problems == 0):
                final_menu = "START"
            elif (count_problems != 0):
                concat_menu = "START"
            input_menu = "QUIT"
    input_menu_complete = molidutils.check_complete(molIDConversion_list)

#########################################
# New menu to finalize without printing #
# concatenation menu                    #
#########################################

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



########################################
# INTERACTIVE MOLID CONCATENATION MENU #
########################################


if (concat_menu == "START"):
    # Prepare for concatenation step
    # We now have to take the information contained in the MolIDConversion objects
    # in molIDConversion_list to update the MolID objects in master_molID_class_list
    # We then need to mine these updated MolID objects to figure out which ones
    # contain concatenated chains. These will be presented to the user in another
    # interactive menu section where they can update the planned conversion on
    # a file by file basis

    master_molID_class_list = molidutils.update_masterlist(master_molID_class_list, molIDConversion_list)

    concat_menu = ""
    concat_menu_complete = ""

    while(concat_menu != "QUIT"):

        count_problems = molidutils.problem_counter(master_molID_class_list)
        if (count_problems == 0):
            concat_menu_complete = "1"

        if (concat_menu_complete == "1"):
            print("""Congratulations! You have successfully constructed your
                     conversion templates.You can proceed to the next section
                     by selection option 6 or, continue to edit your conversion
                     template through this menu
                  """)
        print("""CIFcurator Concatenations Menu
                 -------------------------------
                 This menu appeared because you have assigned the same chain name to two (or more) entities.
                 Note that this will concatenate the entities. So you need to either re-assign chain names,
                 or ACCEPT concatenations.

                 Note: All proposed concatenations must be accepted (by running option 4 or 5) before the curation can
                 be completed.

                 Before you do anything, we suggest to choose option 2, so you know which concatenations have not
                 been accepted. It will also give you the proper format of the input for option 3.

                 If you are sure that all the concatenations are correct. Option 5 will accept all of them. They 
                 will be printed to screen as they are being accepted. 

                 Select one of the following options to proceed:
                 1) Show all conversions
                 2) Show only unaccepted concatenations
                 3) Search and modify destination chainIDs of proposed concatenations
                 4) Accept proposed concatenation one by one
                    (Repeat this step until finalizing option appears)
                 5) Accept ALL (BE CAREFUL, make sure you agree with all concatenations)
              """)
        if (concat_menu_complete == "1"):
            print("    6) Finalize Curation")

        concat_menu = input('Option Number: ')

        if (concat_menu == "1"):
            molidutils.show_full_conversion(master_molID_class_list, step='concatenation')
        elif (concat_menu == "2"):
            molidutils.show_unassigned_conversion(master_molID_class_list, step='concatenation')
        elif (concat_menu == "3"):
            master_molID_class_list = molidutils.edit_concatenation_interface(master_molID_class_list, action='try')[0]
        elif (concat_menu == "4"):
            unassigned_MolID=molidutils.return_unassigned_conversion(master_molID_class_list, step='concatenation')[0]
            print("This is the concatenation you need to accept:")
            new_order=None
            master_molID_class_list = molidutils.list_accept_concatenations(master_molID_class_list, unassigned_MolID, new_order=new_order, action='accept')[0]
            # Note for tomorrow: here we need to create a new function in molidutils, so we can go over all concatenations!
        elif (concat_menu == "5"):
            allnewchains=molidutils.return_unassigned_conversion(master_molID_class_list, step='concatenation')
            for newchain in allnewchains:
                new_order=None
                master_molID_class_list = molidutils.list_accept_concatenations_auto(master_molID_class_list, newchain, new_order=new_order, action='accept')[0]
            count_problems = molidutils.problem_counter(master_molID_class_list)
            print(count_problems)
        elif (concat_menu == "6"):
            print("Finalizing Curation ...")
            molidutils.masterlist_to_pdb(filelist, master_molID_class_list, target_dir=target_dir)
            concat_menu = "QUIT"
