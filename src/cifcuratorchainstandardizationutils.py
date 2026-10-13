from __future__ import print_function
from __future__ import division
import sys
import argparse
import re
import numpy as np
import csv
import os
import copy
import time
from Bio.PDB.MMCIFParser import FastMMCIFParser
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from CIFcurator.alignmentutils import *
from CIFcurator.listutils import *
from matching.games import HospitalResident
import json
sys.setrecursionlimit(10000)

####################
# INITIALIZE STEPS #
####################

def pdb_to_structurelists(filelist):
    """
    Iterates through a list CIF(s) and retrieves structure IDs, chain IDs, and maps chain IDs to their sequences.

    Parameters:
    ----------
    filelist : str
    	list of file paths for all '.cif' files in specified directory

    Returns:
    ------
    Structure_Sequences : list of dict
        Contains dictionaries where each dictionary maps chain IDs to their
        sequences for each structure.
    structid_list : list of str
    	List of unique structure identifiers for each CIF. Format is 'input directory / CIF'
    chid_list : list of str
    	Contains all the chain IDs from CIF(s)
    """
    Structure_Sequences = []
    chid_list = []
    structid_list = []
    N = 0
    for my_file in filelist:
        N += 1
        print("Reading:" + ' ' + my_file + "  (" + str(N) + " of " + str(len(filelist)) + ")")
        struct = FastMMCIFParser(QUIET=1).get_structure(str(my_file), my_file)
        structid_list.append(struct.get_id())
        chid_seq_map = {}
        for chain in struct[0]:
            seq = ""
            for residue in chain:
                seq += ResnConvert(residue.get_resname())
            seq = re.sub('X', '', seq)
            if (len(seq) > 4):
                chid_seq_map[chain.get_id()] = seq
                chid_list.append(chain.get_id())
        Structure_Sequences.append(chid_seq_map)
    chid_set = set(chid_list)
    chid_list = sorted(list(chid_set))
    return Structure_Sequences, structid_list, chid_list


#########################################
# INTERACTIVE STANDARDIZATION FUNCTIONS #
#########################################

def select_standard_seq_from_reference(Structure_Sequences, Standard_Sequences, structid_list, input_menu_check_1):
    """
    User interface for selecting and managing standard sequences from input structures.

    Users can view a list of input structures, select a structure, remove chains from the
    standard sequences, inspect chains, or return to the main menu.

    Parameters:
    -----------
    Structure_Sequences : list of dict
        Contains dictionaries where each dictionary maps chain IDs to their sequences for each structure.
    Standard_Sequences : dict
        Empty dictionary
    structid_list : list of str
        List of unique structure identifiers for each CIF. Format is 'input directory / CIF'
    input_menu_check_1 : str
        Option chosen by user which opens the submenu

    Returns:
    --------
    Standard_Sequences : dict
        The updated dictionary of standard sequences where keys are chain IDs and values are
        their corresponding sequences from the selected CIF. Updated dictionary reflects any
        changes made through user interactions.
    input_menu_check_1 : str
        Updated string representing the state of the main menu, set to '1' to indicate a state change.
    """
    input_submenu = ""
    input_submenu_check_1 = ""
    while(input_submenu != "QUIT"):
        print("    Select Standard Sequences from input structure",
              "    1) Show list of input structures",
              "    2) Select input structure",
              sep="\n")
        if (input_submenu_check_1 == "1"):
            print("    3) Remove chains from Standard Sequences",
                  "    4) Inspect chains in Standard Sequences",
                  "    5) Return to main menu",
                  sep="\n")
        input_submenu = input('Option Number: ')
        if (input_submenu == "1"):
            show_list(structid_list)
        elif (input_submenu == "2"):
            Standard_Sequences, input_menu_check_1, input_submenu_check_1 = select_input_structure(Structure_Sequences,
                                                                                                   structid_list,
                                                                                                   input_menu_check_1,
                                                                                                   input_submenu_check_1)

        elif (input_submenu == "3" and input_submenu_check_1 == "1"):
            Standard_Sequences = remove_chains_from_standard(Standard_Sequences)
        elif (input_submenu == "4" and input_submenu_check_1 == "1"):
            inspect_chains_in_standard(Standard_Sequences)
        elif (input_submenu == "5" and input_submenu_check_1 == "1"):
            input_submenu = "QUIT"
    return Standard_Sequences, input_menu_check_1

def create_standard_seq_from_consensus(Structure_Sequences, Standard_Sequences, chid_list, input_menu_check_1):
    """
    User interface to manage the generation and update of standard sequences for a list of chain IDs
    based on input structures.

    Parameters:
    -----------
    Structure_Sequences : list of dict
        Contains dictionaries where each dictionary maps chain IDs to their sequences for each structure.
    Standard_Sequences : dict
        Empty dictionary
    chid_list : list of str
       Contains all the chain IDs from CIF(s)
    input_menu_check_1 : str
        Option chosen by user which opens the submenu

    Returns:
    --------
    Standard_Sequences : dict
        A dictionary updated with the chain ID as the key and the sequence with the highest consensus
        score as the value.
    input_menu_check_1 : str
        Updated string representing the state of the main menu, set to '1' to indicate a state change.
    """
    input_submenu = ""
    while(input_submenu != "QUIT"):
        print("    Generate Standard Sequences based on all the input structures.",
              "    Type QUIT to return to the main menu.",
              "    1) Show list of chain IDs for Standard Sequences",
              "    2) Enter chain IDs to remove from list",
              "    3) Input file with list of chain IDs to remove",
              "    4) Generate Standard Sequences based on all input structures",
              "    5) Load previously generated standard sequences",
              sep="\n")
        input_submenu = input('Option Number: ')
        if (input_submenu == "1"):
            show_list(chid_list)
        elif (input_submenu == "2"):
            chid_list = remove_user_defined_chain_from_list(chid_list)
        elif (input_submenu == "3"):
            chid_list = remove_file_defined_chain_from_list(chid_list)
        elif (input_submenu == "4"):
            for chid in chid_list:
                print(chid)
                Standard_Sequences = assign_standard_from_consensus(Structure_Sequences, Standard_Sequences, chid)
                print(Standard_Sequences)
            input_submenu = "QUIT"
            input_menu_check_1 = "1"
        elif (input_submenu == "5"):
            Std_Seq_file = input("File with dictionary of the standard sequences: ")
            with open(Std_Seq_file) as f:
                data = f.read()
            print(type(data))
            Standard_Sequences = json.loads(data)
            print(type(Standard_Sequences))

            input_submenu = "QUIT"
            input_menu_check_1 = "1"
    return Standard_Sequences, input_menu_check_1

def review_standard_seq(Structure_Sequences, Standard_Sequences):
    """
    User interface for reviewing and editing standard sequences.

    Users can view chain IDs in the standard sequences, inspect or edit sequences for a
    specific chain ID, and analyze the consensus of sequences for a given chain ID across
    multiple structures

    Parameters:
    -----------
    Structure_Sequences : list of dict
        Contains dictionaries where each dictionary maps chain IDs to their
        sequences for each structure.
    Standard_Sequences : dict
        dictionary where each key is a chain ID and each value is the sequence associated
        with that chain ID, reflecting updates made in previous menu options.

    Returns:
    --------
    None
    """
    input_submenu = ""
    while(input_submenu != "QUIT"):
        print("    Review Standard Sequences. Type QUIT to return to the main menu.",
              "    1) Show list of chain IDs in Standard Sequences",
              "    2) Enter chain ID of Standard Sequence to inspect/edit",
              "    3) Enter chain ID and inspect consensus of all matching chains",
              sep="\n")
        input_submenu = input('Option Number: ')
        if (input_submenu == "1"):
            show_list(Standard_Sequences)
        elif (input_submenu == "2"):
            input_submenu = input('Chain ID: ')
            if input_submenu in Standard_Sequences:
                print(Standard_Sequences[input_submenu])
        elif (input_submenu == "3"):
            chid = input('Chain ID: ')
            this_chainsseq_list, this_chainsseq_score = get_this_chainsseq_list(Structure_Sequences, chid, verbose=True)

def align_to_std_seq_and_save_to_disk(Structure_Sequences, Standard_Sequences, structid_list, filelist, target_dir):
    """
    User interface for performing pairwise alignments of sequences in input structures against standard sequences
    and saves the results.

    Allows users to ignore specific chain IDs during alignment. It manages alignment results, renaming chains if
    necessary, and saving new files to disk.

    Parameters:
    -----------
    Structure_Sequences : list of dict
        Contains dictionaries where each dictionary maps chain IDs to their
        sequences for each structure.
    Standard_Sequences : dict
        Dictionary where each key is a chain ID and each value is the sequence associated
        with that chain ID, reflecting updates made in previous menu options.
    structid_list : list of str
    	List of unique structure identifiers for each CIF. Format is 'input directory / CIF'
    filelist : str
    	list of file paths for all '.cif' files in specified directory
    target_dir : str
        Directory where the new files will be saved

    Returns:
    --------
    None
    """
    ignore_chid = []
    input_submenu = ""
    while(input_submenu != "QUIT"):
        print("    Perform pairwise alignments against Standard Sequences. Type QUIT to return to the main menu.",
              "    1) Show list of structure chain IDs to ignore when pairwise aligning to the Standard Sequences",
              "    2) Enter chain IDs to add to ignore list",
              "    3) Input file with list of chain IDs to add to ignore list",
              "    4) Perform pairwise alignments against Standard Sequences and rename chains (if needed)",
              sep="\n")
        input_submenu = input('Option Number: ')
        if (input_submenu == "1"):
            for chid in sorted(ignore_chid):
                print(chid)
        elif (input_submenu == "2"):
            print("    Enter chain ID to add to ignore list. When complete, enter DONE.")
            while(input_submenu != "DONE"):
                input_submenu = input('Chain ID: ')
                if (input_submenu != "DONE"):
                    ignore_chid.append(input_submenu)
        elif (input_submenu == "3"):
            input_submenu = input('File: ')
            if (os.path.isfile(input_submenu) == True):
                with open(input_submenu) as my_file:
                    for line in my_file.readlines():
                        ignore_chid.append(line.strip())
            else:
                print("File does not exist.")
        elif (input_submenu == "4"):
            counter=0

            for chid_seq_map in Structure_Sequences:
                ChainReassignmentMapping_List = []
                ChainReassignmentScores_List = []

                filelist2=[]
                filelist2.append(filelist[counter])

                print('I am starting to work on:')
                print(filelist2)
                print("this is chid_seq_map and length")
                print(chid_seq_map)
                chid_seq_map_len = len(chid_seq_map)
                print(chid_seq_map_len)
                std_chid_list = []
                test_list_list= {}
                for std_chid in Standard_Sequences:
                    if std_chid not in ignore_chid: # Standard_Sequences is the dictionary with the {chain IDs:sequences} from the reference structure
                        print("Now I am working on this chain:")
                        print(std_chid)
                        std_chid_list.append(std_chid)

                        structchid_score_map = {}
                        if std_chid in chid_seq_map:
                            aligned_seq = AlignSequences([Standard_Sequences[std_chid], chid_seq_map[std_chid]])
                            score = ScoreSequenceAlignment(aligned_seq[0], aligned_seq[1])
                            # If score is perfect, don't bother aligning all chains
                            if (score >= 0.85):
                                structchid_score_map[std_chid] = score
                                for struct_chid in chid_seq_map:
                                    if struct_chid not in ignore_chid:
                                        if (struct_chid != std_chid):
                                            structchid_score_map[struct_chid] = 0
                        # If score is not perfect, align all chains
                            else:
                                for struct_chid in chid_seq_map:
                                    if struct_chid not in ignore_chid:
                                        aligned_seq = AlignSequences([Standard_Sequences[std_chid], chid_seq_map[struct_chid]])
                                        score = ScoreSequenceAlignment(aligned_seq[0], aligned_seq[1])
                                        structchid_score_map[struct_chid] = score
                    # Chid of standard does not match any structure chid, so perform all alignments
                        else:
                            for struct_chid in chid_seq_map:
                                if struct_chid not in ignore_chid:
                                    aligned_seq = AlignSequences([Standard_Sequences[std_chid], chid_seq_map[struct_chid]])
                                    score = ScoreSequenceAlignment(aligned_seq[0], aligned_seq[1])
                                    structchid_score_map[struct_chid] = score


                    #if std_chid not in ignore_chid: # we don't want to have any of the ignored chains 0_0
                        test_list_list[std_chid] = structchid_score_map

                        print("Chains already completed:")
                        print(std_chid_list)

                ######
                #  In this section we re-structure the data we collected in previous step. In previous step we filled in a
                #  "matrix" of chainids and scores. In this section we first create the 'transpose', and then we sort both
                # dictionaries of chainids by their scores.
                # This step is necessary to later use the matching algorithm to assign the new chain ids.
                #
                ######
                test_list_list_T = {}
                for tll_key in test_list_list.keys():
                    for tll_sub_key in test_list_list[tll_key].keys():
                        if tll_sub_key not in test_list_list_T:
                            test_list_list_T[tll_sub_key] = {}
                        test_list_list_T[tll_sub_key][tll_key] = test_list_list[tll_key][tll_sub_key]


                test_list_list_sortkey2 = {}
                for tll_key in test_list_list.keys():
                    test_list_list_sortkey2[tll_key] = sorted(test_list_list[tll_key], key=lambda x:test_list_list[tll_key][x], reverse=True)


                test_list_list_T_sortkey = {}
                for tll_key in test_list_list_T.keys():
                    test_list_list_T_sortkey[tll_key] = sorted(test_list_list_T[tll_key], key=lambda x:test_list_list_T[tll_key][x], reverse=True)

                ########
                # The following section uses the matching python library to assign the chain IDs. We use the HospitalResident algorithm
                # We create a capacities variable, that assigns only one chain per hospital.
                ########

                capacities={}
                for h in list(test_list_list_sortkey2.keys()):
                    capacities[h] = 1

                game = HospitalResident.create_from_dictionaries(test_list_list_sortkey2, test_list_list_T_sortkey, capacities)
                solved_game = game.solve()

                output = {}
                for k in solved_game.keys():
                    output[str(k)] = str(solved_game[k][0])

                output_scores = {}
                for k in output.keys():
                    output_scores[str(k)] = str(test_list_list_T[k][output[k]])


                print('Just finished with this structure:')
                print(filelist2)
                print('These are the results:')
                print(output)
                print(output_scores)

                ChainReassignmentMapping_List.append(output) # FAPA NEW SCORES
                ChainReassignmentScores_List.append(output_scores) # FAPA NEW SCORES

                ##
                # Here is where we start writing the new structures
                ##

                reassignedmaps_to_pdb(filelist2, ChainReassignmentMapping_List, filelist2, target_dir=target_dir)
                reassignedmaps_to_log(ChainReassignmentMapping_List, ChainReassignmentScores_List, filelist2, target_dir=target_dir)

                new_pdb_out=target_dir+"/"+(filelist2[0]).split('/')[-1]
                print("new_pdb_out")
                print(new_pdb_out)

                while not os.path.exists(new_pdb_out):
                    time.sleep(1) #FAPA, WAITING LESS TIME
                    print("waiting...")

                while not os.path.getsize(new_pdb_out) > 0:
                    time.sleep(1) #FAPA, WAITING LESS TIME
                    print("waiting even more...")

                counter+=1

                input_submenu = "QUIT"

 # Not called anywhere
def align_to_standard_seq(Structure_Sequences, Standard_Sequences, structid_list):
    """
    User interface for performing pairwise alignments against Standard Sequences.

    Parameters:
    -----------
    Structure_Sequences : list of dict
        Contains dictionaries where each dictionary maps chain IDs to their
        sequences for each structure.
    Standard_Sequences : dict
        dictionary where each key is a chain ID and each value is the sequence associated
        with that chain ID, reflecting updates made in previous menu options.
    structid_list : list of str
    	List of unique structure identifiers for each CIF. Format is 'input directory / CIF'

    Returns:
    --------
    ChainReassignmentMapping_List : list of dict
        A list of dictionaries where each dictionary maps original chain IDs to reassigned chain IDs.
    ChainReassignmentScores_List : list of dict
        A list of dictionaries where each dictionary contains scores for the chain reassignment mapping.
    input_menu_check_2 : str
        Status indicator for menu navigation, typically set to "1" to continue.
    """
    ignore_chid = []
    input_submenu = ""
    while(input_submenu != "QUIT"):
        print("    Perform pairwise alignments against Standard Sequences. Type QUIT to return to the main menu.",
              "    1) Show list of structure chain IDs to ignore when pairwise aligning to the Standard Sequences",
              "    2) Enter chain IDs to add to ignore list",
              "    3) Input file with list of chain IDs to add to ignore list",
              "    4) Perform pairwise alignments against Standard Sequences and create conversion template (MAY NEED USER INPUT)",
              "    5) Perform pairwise alignments against Standard Sequences and AUTOMATICALLY create conversion template",
              "    6) TESTING",
              sep="\n")
        input_submenu = input('Option Number: ')
        if (input_submenu == "1"):
            for chid in sorted(ignore_chid):
                print(chid)
        elif (input_submenu == "2"):
            print("    Enter chain ID to add to ignore list. When complete, enter DONE.")
            while(input_submenu != "DONE"):
                input_submenu = input('Chain ID: ')
                if (input_submenu != "DONE"):
                    ignore_chid.append(input_submenu)
        elif (input_submenu == "3"):
            input_submenu = input('File: ')
            if (os.path.isfile(input_submenu) == True):
                with open(input_submenu) as my_file:
                    for line in my_file.readlines():
                        ignore_chid.append(line.strip())
            else:
                print("File does not exist.")
        elif (input_submenu == "4"):
            ChainReassignmentMapping_List = []
            ChainReassignmentScores_List = []
            for chid_seq_map in Structure_Sequences:
                Struct_ChainReassignmentMap = {}
                Struct_ChainReassignmentScores = {}
                std_chid_list = []
                chidlist_list = []
                scorelist_list = []
                for std_chid in Standard_Sequences:
                    std_chid_list.append(std_chid)
                    structchid_score_map = {}
                    # First, check to see if chid is already correctly assigned
                    if std_chid in chid_seq_map:
                        aligned_seq = AlignSequences([Standard_Sequences[std_chid], chid_seq_map[std_chid]])
                        score = ScoreSequenceAlignment(aligned_seq[0], aligned_seq[1])
                        # If score is perfect, don't bother aligning all chains
                        if (score == 1):
                            structchid_score_map[std_chid] = score
                            for struct_chid in chid_seq_map:
                                if (struct_chid != std_chid):
                                    structchid_score_map[struct_chid] = 0
                        # If score is not perfect, align all chains
                        else:
                            for struct_chid in chid_seq_map:
                                if struct_chid not in ignore_chid:
                                    aligned_seq = AlignSequences([Standard_Sequences[std_chid], chid_seq_map[struct_chid]])
                                    score = ScoreSequenceAlignment(aligned_seq[0], aligned_seq[1])
                                    structchid_score_map[struct_chid] = score
                    # Chid of standard does not match any structure chid, so perform all alignments
                    else:
                        for struct_chid in chid_seq_map:
                            if struct_chid not in ignore_chid:
                                aligned_seq = AlignSequences([Standard_Sequences[std_chid], chid_seq_map[struct_chid]])
                                score = ScoreSequenceAlignment(aligned_seq[0], aligned_seq[1])
                                structchid_score_map[struct_chid] = score
                    scorelist_list.append(list(reversed(sorted(structchid_score_map.values()))))
                    chidlist_list.append(list(reversed(sorted(structchid_score_map, key=structchid_score_map.get))))
                # Relevant information to create chain reassingment map has now been gathered from each structure via alignment
                # std_chid_list, chidlist_list and scorelist_list have the same indexes
                # At a given location i, std_chid_list[i] is the standard_sequence chid,
                # chidlist_list[i] is the list of chain ids in the given structure aligned to that standard sequence
                # and, scorelist_list[i] is the list of corresponding scores of that alignment
                # Note that there is no reason to keep all scores of all chids instead of just the top
                # scoring one, however, it's done because it seems likely that information may
                # be valuable and worth incorporating later on
                # Construct the map, making sure to avoid conflicts by comparing scores
                # These maps are old to new
                # CHANGES MADE HERE TO MAKE SURE CONFLICTS DO NOT ARRISE BETWEEN UNASSIGNED CHAINIDs and
                for i in range(len(std_chid_list)):
                    print("This is the dictionary at beginning of loop:")
                    print(Struct_ChainReassignmentMap)
                    print("This chidlist_list[i][0]:")
                    print(chidlist_list[i][0])
                    print("This is std_chid_list[i]:")
                    print(std_chid_list[i])
                    if chidlist_list[i][0] not in Struct_ChainReassignmentMap:
                        Struct_ChainReassignmentMap[chidlist_list[i][0]] = std_chid_list[i]
                        Struct_ChainReassignmentScores[chidlist_list[i][0]] = scorelist_list[i][0]
                    else:
                        # I < i
                        I = std_chid_list.index(Struct_ChainReassignmentMap[chidlist_list[i][0]])
                        # Compare score of original and new entry. If new entry has better score, replace
                        if(scorelist_list[i][0] > scorelist_list[I][0]):
                            Struct_ChainReassignmentMap[chidlist_list[i][0]] = std_chid_list[i]
                            Struct_ChainReassignmentScores[chidlist_list[i][0]] = scorelist_list[i][0]
                            if scorelist_list[I][0] in Struct_ChainReassignmentScores:
                                del Struct_ChainReassignmentScores[chidlist_list[I][0]]
                                del Struct_ChainReassignmentMap[chidlist_list[I][0]]

                    print("This is the dictionary at end of loop:")
                    print(Struct_ChainReassignmentMap)

                # Now for final check, make sure no conflicts by checking unused chid
                # Loop over all chid in the structure
                unused_chid = []
                free_chid = [] #FAPA
                for chid in chid_seq_map:
                    if chid not in Struct_ChainReassignmentMap.keys():
                        unused_chid.append(chid)
                        # chid is unused . . . will not change. If it is in destination list, then there will be conflict
                    if chid not in Struct_ChainReassignmentMap.values():
                        free_chid.append(chid)
                for chid in unused_chid:
                    used_dest_chid = list(Struct_ChainReassignmentMap.values()) #FAPA
                    if chid in used_dest_chid:
                        # chid is unused and a destination so a conflict exists and needs to be given destination
                        II = Structure_Sequences.index(chid_seq_map)
                        print("Conflict in structure " + structid_list[II])
                        for old_chid in Struct_ChainReassignmentMap:
                            if (Struct_ChainReassignmentMap[old_chid] == chid):
                                valid_chid = False #LIV
                                while valid_chid == False: #LIV
                                    print("Original chain " + str(old_chid) + " is being reassigned to new chain " + str(chid))
                                    print("This creates a conflict because original chain " + str(chid) + " has not been reassigned.")
                                    print("Where should original chain " + str(chid) + " be reassigned?")
                                    print("You can give a new id or choose from:") #FAPA 011823
                                    print(free_chid)
                                    user_input_submenu = input('New chain ID: ') #FAPA
                                    if (user_input_submenu not in unused_chid) and (user_input_submenu not in Struct_ChainReassignmentMap.values()): #FAPA
                                        # input is a valid destination
                                        valid_chid = True
                                    else:
                                        print("Not a valid destination chain ID")
                        if (len(input_submenu) <= 2):
                            Struct_ChainReassignmentMap[chid] = user_input_submenu
                            Struct_ChainReassignmentScores[chid] = 0
                            used_dest_chid.append(user_input_submenu)
                        else:
                            Struct_ChainReassignmentMap[chid] = user_input_submenu[0:2]
                            Struct_ChainReassignmentScores[chid] = 0
                            used_dest_chid.append(user_input_submenu[0:2])
                ChainReassignmentMapping_List.append(Struct_ChainReassignmentMap)
                ChainReassignmentScores_List.append(Struct_ChainReassignmentScores)
                input_menu_check_2 = "1"
                input_submenu = "QUIT"

#FAPA ADDS A NEW Option

        elif (input_submenu == "5"):
            ChainReassignmentMapping_List = []
            ChainReassignmentScores_List = []
            for chid_seq_map in Structure_Sequences:
                Struct_ChainReassignmentMap = {}
                Struct_ChainReassignmentScores = {}
                std_chid_list = []
                chidlist_list = []
                scorelist_list = []
                for std_chid in Standard_Sequences:
                    std_chid_list.append(std_chid)
                    structchid_score_map = {}
                    # First, check to see if chid is already correctly assigned
                    if std_chid in chid_seq_map:
                        aligned_seq = AlignSequences([Standard_Sequences[std_chid], chid_seq_map[std_chid]])
                        score = ScoreSequenceAlignment(aligned_seq[0], aligned_seq[1])
                        # If score is perfect, don't bother aligning all chains
                        if (score == 1):
                            structchid_score_map[std_chid] = score
                            for struct_chid in chid_seq_map:
                                if (struct_chid != std_chid):
                                    structchid_score_map[struct_chid] = 0
                        # If score is not perfect, align all chains
                        else:
                            for struct_chid in chid_seq_map:
                                if struct_chid not in ignore_chid:
                                    aligned_seq = AlignSequences([Standard_Sequences[std_chid], chid_seq_map[struct_chid]])
                                    score = ScoreSequenceAlignment(aligned_seq[0], aligned_seq[1])
                                    structchid_score_map[struct_chid] = score
                    # Chid of standard does not match any structure chid, so perform all alignments
                    else:
                        for struct_chid in chid_seq_map:
                            if struct_chid not in ignore_chid:
                                aligned_seq = AlignSequences([Standard_Sequences[std_chid], chid_seq_map[struct_chid]])
                                score = ScoreSequenceAlignment(aligned_seq[0], aligned_seq[1])
                                structchid_score_map[struct_chid] = score
                    scorelist_list.append(list(reversed(sorted(structchid_score_map.values()))))
                    chidlist_list.append(list(reversed(sorted(structchid_score_map, key=structchid_score_map.get))))
                # Relevant information to create chain reassingment map has now been gathered from each structure via alignment
                # std_chid_list, chidlist_list and scorelist_list have the same indexes
                # At a given location i, std_chid_list[i] is the standard_sequence chid,
                # chidlist_list[i] is the list of chain ids in the given structure aligned to that standard sequence
                # and, scorelist_list[i] is the list of corresponding scores of that alignment
                # Note that there is no reason to keep all scores of all chids instead of just the top
                # scoring one, however, it's done because it seems likely that information may
                # be valuable and worth incorporating later on
                # Construct the map, making sure to avoid conflicts by comparing scores
                # These maps are old to new
                # CHANGES MADE HERE TO MAKE SURE CONFLICTS DO NOT ARRISE BETWEEN UNASSIGNED CHAINIDs and
                for i in range(len(std_chid_list)):
                    print("This is the dictionary at beginning of loop:")
                    print(Struct_ChainReassignmentMap)
                    print("This chidlist_list[i][0]:")
                    print(chidlist_list[i][0])
                    print("This is std_chid_list[i]:")
                    print(std_chid_list[i])
                    if chidlist_list[i][0] not in Struct_ChainReassignmentMap:
                        Struct_ChainReassignmentMap[chidlist_list[i][0]] = std_chid_list[i]
                        Struct_ChainReassignmentScores[chidlist_list[i][0]] = scorelist_list[i][0]
                    else:
                        # I < i
                        I = std_chid_list.index(Struct_ChainReassignmentMap[chidlist_list[i][0]])
                        # Compare score of original and new entry. If new entry has better score, replace
                        if(scorelist_list[i][0] > scorelist_list[I][0]):
                            Struct_ChainReassignmentMap[chidlist_list[i][0]] = std_chid_list[i]
                            Struct_ChainReassignmentScores[chidlist_list[i][0]] = scorelist_list[i][0]
                            if scorelist_list[I][0] in Struct_ChainReassignmentScores:
                                del Struct_ChainReassignmentScores[chidlist_list[I][0]]
                                del Struct_ChainReassignmentMap[chidlist_list[I][0]]

                    print("This is the dictionary at end of loop:")
                    print(Struct_ChainReassignmentMap)

                # Now for final check, make sure no conflicts by checking unused chid
                # Loop over all chid in the structure
                unused_chid = []
                free_chid=[] #FAPA
                counter=0 #FAPA
                for chid in chid_seq_map:
                    if chid not in Struct_ChainReassignmentMap.keys():
                        unused_chid.append(chid)
                        # chid is unused . . . will not change. If it is in destination list, then there will be conflict
                    if chid not in Struct_ChainReassignmentMap.values():
                        free_chid.append(chid)
                for chid in unused_chid:
                    used_dest_chid = list(Struct_ChainReassignmentMap.values()) #FAPA
                    if chid in used_dest_chid:
                        # chid is unused and a destination so a conflict exists and needs to be given destination
                        II = Structure_Sequences.index(chid_seq_map)
                        print("Conflict in structure " + structid_list[II])
                        for old_chid in Struct_ChainReassignmentMap:
                            if (Struct_ChainReassignmentMap[old_chid] == chid):
                                user_input_submenu = free_chid[counter]
                                counter+=1

                        Struct_ChainReassignmentMap[chid] = user_input_submenu
                        Struct_ChainReassignmentScores[chid] = 0
                        used_dest_chid.append(user_input_submenu)

                ChainReassignmentMapping_List.append(Struct_ChainReassignmentMap)
                ChainReassignmentScores_List.append(Struct_ChainReassignmentScores)
                input_menu_check_2 = "1"
                input_submenu = "QUIT"

#FAPA FINISHES ADDING NEW OPTION

#FAPA ADDS OPTION 6

        elif (input_submenu == "6"):
            ChainReassignmentMapping_List = []
            ChainReassignmentScores_List = []
            for chid_seq_map in Structure_Sequences:
                print('I am starting to work on:')
                print(structid_list[Structure_Sequences.index(chid_seq_map)])
                print("this is chid_seq_map and length")
                print(chid_seq_map)
                chid_seq_map_len=len(chid_seq_map)
                print(chid_seq_map_len)
                std_chid_list = []
                test_list_list= {}
                for std_chid in Standard_Sequences:
                    if std_chid not in ignore_chid: # Standard_Sequences is the dictionary with the {chain IDs:sequences} from the reference structure
                        print("this is std_chid:")
                        print(std_chid)

                        print("IGNORE IGNORE IGNORE:")
                        print(ignore_chid)

                        std_chid_list.append(std_chid)

                        structchid_score_map = {}
                        # First, check to see if chid is already correctly assigned
                        if std_chid in chid_seq_map:
                        #if std_chid not in ignore_chid: # FAPA
                            aligned_seq = AlignSequences([Standard_Sequences[std_chid], chid_seq_map[std_chid]])
                            score = ScoreSequenceAlignment(aligned_seq[0], aligned_seq[1])
                            # If score is perfect, don't bother aligning all chains
                            if (score >= 0.85): # No need to be so strict, think of possible cases when we would need identical score?
                                structchid_score_map[std_chid] = score
                                for struct_chid in chid_seq_map:
                                    if struct_chid not in ignore_chid:
                                        if (struct_chid != std_chid):
                                            structchid_score_map[struct_chid] = 0
                        # If score is not perfect, align all chains
                            else:
                                for struct_chid in chid_seq_map:
                                    if struct_chid not in ignore_chid:
                                        aligned_seq = AlignSequences([Standard_Sequences[std_chid], chid_seq_map[struct_chid]])
                                        score = ScoreSequenceAlignment(aligned_seq[0], aligned_seq[1])
                                        structchid_score_map[struct_chid] = score
                    # Chid of standard does not match any structure chid, so perform all alignments
                        else:
                            for struct_chid in chid_seq_map:
                                if struct_chid not in ignore_chid:
                                    aligned_seq = AlignSequences([Standard_Sequences[std_chid], chid_seq_map[struct_chid]])
                                    score = ScoreSequenceAlignment(aligned_seq[0], aligned_seq[1])
                                    structchid_score_map[struct_chid] = score


                    #if std_chid not in ignore_chid: # we don't want to have any of the ignored chains 0_0
                        test_list_list[std_chid]= structchid_score_map

                        print("test_list_list:")
                        print(test_list_list)
                        print("std_chid_list:")
                        print(std_chid_list)



                ######
                #  In this section we re-structure the data we collected in previous step. In previous step we filled in a
                #  "matrix" of chainids and scores. In this section we first create the 'transpose', and then we sort both
                # dictionaries of chainids by their scores.
                # This step is necessary to later use the matching algorithm to assign the new chain ids.
                #
                ######
                test_list_list_T = {}
                for tll_key in test_list_list.keys():
                    for tll_sub_key in test_list_list[tll_key].keys():
                        if tll_sub_key not in test_list_list_T:
                            test_list_list_T[tll_sub_key] = {}
                        test_list_list_T[tll_sub_key][tll_key] = test_list_list[tll_key][tll_sub_key]


                test_list_list_sortkey2 = {}
                for tll_key in test_list_list.keys():
                    test_list_list_sortkey2[tll_key] = sorted(test_list_list[tll_key], key=lambda x:test_list_list[tll_key][x], reverse=True)
                print('test_list_list_sortkey2')
                print(test_list_list_sortkey2)


                test_list_list_T_sortkey = {}
                for tll_key in test_list_list_T.keys():
                    test_list_list_T_sortkey[tll_key] = sorted(test_list_list_T[tll_key], key=lambda x:test_list_list_T[tll_key][x], reverse=True)
                print('test_list_list_T_sortkey')
                print(test_list_list_T_sortkey)

                print("std_chid_list")
                print(len(std_chid_list))

                print(chid_seq_map_len)

                print("len test_list_list_sortkey2.keys")
                print(list(test_list_list_sortkey2.keys()))

                ########
                # The following section uses the matching python library to assign the chain IDs. We use the HospitalResident algorithm
                # We create a capacities variable, that assigns only one chain per hospital.
                ########

                capacities={}
                for h in list(test_list_list_sortkey2.keys()):
                    capacities[h]=1

                game = HospitalResident.create_from_dictionaries(test_list_list_sortkey2, test_list_list_T_sortkey, capacities)
                solved_game=game.solve()

                output = {}
                for k in solved_game.keys():
                    output[str(k)] = str(solved_game[k][0])

                output_scores = {}
                for k in output.keys():
                    output_scores[str(k)] = str(test_list_list_T[k][output[k]])


                print("Just finished with this structure:")
                print(structid_list[Structure_Sequences.index(chid_seq_map)])
                print('These are the results:')
                print(output)
                print(output_scores)

                ChainReassignmentMapping_List.append(output) # FAPA NEW SCORES
                ChainReassignmentScores_List.append(output_scores) # FAPA NEW SCORES
                input_menu_check_2 = "1"
                input_submenu = "QUIT"

#END OPTION 6

    print(ChainReassignmentMapping_List, ChainReassignmentScores_List)
    return ChainReassignmentMapping_List, ChainReassignmentScores_List, input_menu_check_2


def assign_standard_from_consensus(Structure_Sequences, Standard_Sequences, chid):
    """
    Determine the standard sequence for a given chain ID based on consensus across multiple structures.

    Parameters:
    -----------
    Structure_Sequences : list of dict
        Contains dictionaries where each dictionary maps chain IDs to their sequences for each structure.
    Standard_Sequences : dict
        Empty dictionary
    chid : str
        Chain ID from list containing all chain IDs across CIF(s)

    Returns:
    --------
    Standard_Sequences : dict
        A dictionary updated with the chain ID as the key and the sequence with the highest consensus
        score as the value.
    """
    this_chainsseq_list, this_chainsseq_score = get_this_chainsseq_list(Structure_Sequences, chid)
    # Now you have an ordered by length list of sequences pertaining to particular chID_list
    # You want to grab the longest one with the higest score and save it as this chains standard sequences
    this_score_vec = []
    for seq in this_chainsseq_list:
        this_score_vec.append(this_chainsseq_score[seq])
    this_score_vec = sorted(this_score_vec)
    topscore = this_score_vec[len(this_score_vec)-1]
    for seq in this_chainsseq_list:
        if (topscore == this_chainsseq_score[seq]):
            Standard_Sequences[chid] = seq
            break
    return Standard_Sequences

def get_this_chainsseq_list(Structure_Sequences, chid, verbose=True):
    """
    Extract and analyze sequences associated with a specific chain ID from 'Structure_Sequences'.

    Parameters:
    -----------
    Structure_Sequences : list of dict
        Contains dictionaries where each dictionary maps chain IDs to their sequences for each structure.
    chid : str
        Chain ID from list containing all chain IDs across CIF(s)
    verbose : Bool, optional
        If True, prints additional information about the sequences and their counts. Default is True.
    Returns:
    --------
    this_chainsseq_list : list of str
        A list of sequences associated with the specified chain ID, sorted in descending order
        by sequence length.
    this_chainsseq_score : dict
        A dictionary where sequences are mapped with their counts across the structures,
        adjusted to include counts from shorter sequences that are part of longer ones.
    """
    this_chainsseq_list = []
    # Here, we're going through each structure in Structure_Sequences, which holds chid to seq maps
    for chid_seq_map in Structure_Sequences:
        if chid in chid_seq_map:
            this_chainsseq_list.append(chid_seq_map[chid])
    # Now we have a vector containing all the sequences pertaining to chid we want
    # Below is CRAP!
    # this_chainsseq_list is now sorted by length.
    # Start at longest, see if next longest is containing, if so add
    this_chainsseq_score = {}
    # Initialize scores
    this_chainsseq_set = set(this_chainsseq_list)

    print("this_chansseq_set")
    print(this_chainsseq_set)

    for seq in this_chainsseq_set:
        this_chainsseq_score[seq] = 0
    # For each unique seq, count number of times in list
    for seq in this_chainsseq_list:
        this_chainsseq_score[seq] += 1
    # Reduce list to uniq entries
    this_chainsseq_list = sorted(list(this_chainsseq_set))
    this_chainsseq_list.sort(key=len)
    this_chainsseq_list.reverse()

    print("this is after finding uniques:")
    print(this_chainsseq_list)
    # Determine containedness of a sequence within a longer sequence
    for i in range(len(this_chainsseq_list)):
        for j in range(len(this_chainsseq_list)-i-1):
            if(this_chainsseq_list[i+1+j] in this_chainsseq_list[i]):
                this_chainsseq_score[this_chainsseq_list[i]] += this_chainsseq_score[this_chainsseq_list[i+1+j]]
    if verbose:
        for seq in this_chainsseq_list:
            print("Number of structures: " + str(this_chainsseq_score[seq]))
            print(seq)
    return this_chainsseq_list, this_chainsseq_score

def inspect_chains_in_standard(Standard_Sequences):
    """
    Interactive user interface which allows user to inspect the chain IDs in the standard sequences

    Parameters:
    -----------
    Standard_Sequences : dict
        A dictionary where keys are chain IDs and values are their corresponding sequences
        for the selected CIF file. The dictionary is retrieved based on the user's input for
        the CIF file name.

    Returns:
    --------
    None
    """
    input_submenu = ""
    while(input_submenu != "QUIT"):
        print("    Inspect chains in Standard Sequences. Type QUIT to return to the main menu.",
              "    1) Show list of chain IDs in Standard Sequences",
              "    2) Enter chain IDs inspect",
              sep="\n")
        input_submenu = input('Option Number: ')
        if (input_submenu == "1"):
            for chid in Standard_Sequences:
                print(chid)
        elif (input_submenu == "2"):
            input_submenu = input('Chain ID: ')
            if input_submenu in Standard_Sequences:
                print(Standard_Sequences[input_submenu])
            else:
                print(input_submenu + " is not a chain ID in Standard Sequences")

def remove_chains_from_standard(Standard_Sequences):
    """
    Interactive user interface for removing chain IDs from the standard sequences.

    Parameters:
    -----------
    Standard_Sequences : dict
        A dictionary where keys are chain IDs and values are their corresponding sequences
        for the selected CIF file. The dictionary is retrieved based on the user's input for
        the CIF file name.

    Returns:
    --------
    Standard_Sequences : dict
        Updated dictionary where the specified chain ID(s) has been removed from the
        dictionary.
    """
    input_submenu = ""
    while(input_submenu != "QUIT"):
        print("    Remove chains from Standard Sequences. Type QUIT to return to the main menu.",
              "    1) Show list of chain IDs in Standard Sequences",
              "    2) Enter chain IDs to remove from list",
              "    3) Input file with list of chain IDs to remove",
              sep="\n")
        input_submenu = input('Option Number: ')
        if (input_submenu == "1"):
            show_list(Standard_Sequences)
        elif (input_submenu == "2"):
            Standard_Sequences = remove_user_defined_chain_from_list(Standard_Sequences)
        elif (input_submenu== "3"):
            Standard_Sequences = remove_file_defined_chain_from_list(Standard_Sequences)
    return Standard_Sequences
#Dictionary of the structure ID which contains the chain ID and its associated sequence.

def select_input_structure(Structure_Sequences, structid_list, input_menu_check_1, input_submenu_1_check_1):
    """
    Prompts user to enter name of CIF in format 'input directory/CIF' and selects the corresponding
    structure data.

    Parameters:
    ----------
    Structure_Sequences : dict
    	Contains dictionaries where each dictionary maps chain IDs to their sequences for each structure.
    structid_list : list of str
        List of unique structure identifiers for each CIF. Format is 'input directory / CIF'
    input_menu_check_1 : str
    	Tracks the current state of the main menu.
    input_submenu_1_check_1 : str
    	Tracks the current state of the submenu. It determines whether to show additional menu options or not.

    Returns:
    --------
    Standard_Sequences : dict
        A dictionary where keys are chain IDs and values are their corresponding sequences for the selected CIF file.
        The dictionary is retrieved based on the user's input for the CIF.
    input_menu_check_1 : str
    	Updated string representing the state of the main menu, set to '1' to indicate a state change.
    input_submenu_1_check_1 : str
    	Updated string representing the state of the submenu, set to '1' to indicate a state change.
    """
    input_submenu = ""
    while(input_submenu != "QUIT"):
        input_submenu = input('File name: ')
        if input_submenu in structid_list:
            Standard_Sequences = Structure_Sequences[structid_list.index(input_submenu)]
            input_submenu = "QUIT"
            input_submenu_1_check_1 = str(1)
            input_menu_check_1 = str(1)
        else:
            print("""File name entered was not in list of files.""")
    return Standard_Sequences, input_menu_check_1, input_submenu_1_check_1


#################
# FINALIZE STEP #
#################

def reassignedmaps_to_pdb(filelist, ChainReassignmentMapping_List, structid_list, target_dir=None):
    """
    Reassigns chain IDs in CIF files based on a provided mapping and writes the modified files to a
    target directory.

    Parameters:
    -----------
    filelist : str
        list of file paths for all '.cif' files in specified directory
    ChainReassignmentMapping_List : list of dict
        A list of dictionaries where each dictionary maps original chain IDs to reassigned chain IDs.
    structid_list : list of str
    	List of unique structure identifiers for each CIF. Format is 'input directory / CIF'
    target_dir : str, optional
        Directory where the new files will be saved

    Returns:
    --------
    None
    """
    print("structid_list")
    print(structid_list)

    for my_files in filelist:
        newciffilename=target_dir+'/'+my_files.split('/')[-1]
        with open(my_files) as myfile:
            print("my file name")
            print(myfile.name)
            with open(newciffilename, 'w') as newciffile:
                # Now figure out which file is which template
                I = structid_list.index(myfile.name)
                for line in myfile:
                    if (line[0:4] == "ATOM"):
                            # Chains outside map should not exist but just in case
                        line_split = line.split()
                        if line_split[17] in ChainReassignmentMapping_List[I]:
                            newline = line_split[0] + " " + line_split[1] + " " + line_split[2] + " " + line_split[3] + " " + line_split[4] + " " + line_split[5] + " " + ChainReassignmentMapping_List[I][line_split[17]] + " " + line_split[7] + " " + line_split[8] + " " + line_split[9] + " " + line_split[10] + " " + line_split[11] + " " + line_split[12] + " " + line_split[13] + " " + line_split[14] + " " + line_split[15] + " " + line_split[16] + " " + ChainReassignmentMapping_List[I][line_split[17]] + " " + line_split[18] + " " + line_split[19] + "\n" #FAPA TEST
                            newciffile.write(newline)
                        else:
                            newciffile.write(line)
                    else:
                        newciffile.write(line)

def reassignedmaps_to_log(ChainReassignmentMapping_List, ChainReassignmentScores_List, structid_list, target_dir=None, verbose=True):
    """
    Logs chain reassignment information to a file and optionally prints a summary to the console.

    Parameters:
    -----------
    ChainReassignmentMapping_List : list of dict
        A list of dictionaries where each dictionary maps original chain IDs to reassigned chain IDs.
    ChainReassignmentScores_List : list of dict
        A list of dictionaries where each dictionary contains scores for the chain reassignment mapping.
    structid_list : list of str
    	List of unique structure identifiers for each CIF. Format is 'input directory / CIF'
    target_dir : str, optional
        Directory where the new files will be saved
    verbose : bool, optional
        If True, prints a summary of the chain reassignments and their scores to the console. Default is True.
    Returns:
    --------
    None
    """
    logfilename=target_dir+'/ChainStandardizationRecord.txt'
    with open(logfilename, 'a') as recordfile:
        for I in range(len(structid_list)):
            for chid in ChainReassignmentMapping_List[I]:
                recordfile.write(structid_list[I] + ":" + chid + ":" + ChainReassignmentMapping_List[I][chid] + ":" + str(ChainReassignmentScores_List[I][chid]) + "\n")
    if verbose:
        print("Original Chain ID | New Chain ID | Naive Similarity Score ")
        for i in range(len(ChainReassignmentMapping_List)):
            ChainReassignmentMap = ChainReassignmentMapping_List[i]
            ChainReassignmentScore = ChainReassignmentScores_List[i]
            print(structid_list[i])
            for oldchid in ChainReassignmentMap:
                print(oldchid + " " + ChainReassignmentMap[oldchid] + " " + str(ChainReassignmentScore[oldchid]))
