#!/Users/fatima/anaconda3/envs/PDBCleanV2/bin/python
# coding: utf-8

from __future__ import print_function
from __future__ import division
import sys, glob
from CIFcurator import cifcuratorresiduestandardizationutils as resstd

########################
# READ INPUT ARGUMENTS #
########################
n_arg = len(sys.argv)
if(n_arg<3):
    print('Usage error: {0} <source directory> <target directory>'.format(sys.argv[0]))
    sys.exit()
source_dir=sys.argv[1]
target_dir=sys.argv[2]


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
