from __future__ import print_function
import sys
import re
import numpy as np
import csv
import os
import copy
from Bio.PDB.MMCIF2Dict import MMCIF2Dict


####################
# INITIALIZE STEPS #
####################

def pdb_to_masterlist(filelist):
    """
    Reads a list of CIF(s) and process them to make list of molID classes

    Parameters:
    -----------
    filelist : list of str
    Path to read CIF(s)

    Returns:
    -----------
    master_molID_class_list : list
    A list containing molID class objects for each CIF processed.
    """
    master_molID_class_list = []
    N=0
    for my_files in filelist:
        N += 1
        print("Reading:"+' '+my_files+"  ("+str(N)+" of "+str(len(filelist))+")")
        with open(my_files) as myfile:
            this_molID_class = make_MolID_cif(myfile)
            master_molID_class_list.append(this_molID_class)
    return master_molID_class_list

def uniquelist_to_conversionlist(unique_map):
    """
    Converts a unique mapping dictionary into a list of MolIDConversion objects.

    Parameters:
    -----------
    unique_map : dict
        Keys are entity name strings, & Values are lists of chain IDs associated with each entity name.

    Returns:
    --------
    molIDConversion_list : list
        A list of initialized `MolIDConversion` objects created from the unique mapping.
    """
    molIDConversion_list = []
    for key in unique_map:
        molIDConversion_list.append(initial_MolIDConversion(key, unique_map[key]))
    return molIDConversion_list

def update_masterlist(master_molID_class_list, molIDConversion_list):
    """
    Updates the master list of MolID objects with new chain ID mappings from a list of MolIDConversion objects.

    Function takes the current master list of MolID objects and updates their chain ID mappings using the
    mappings provided in a list of MolIDConversion objects. It also checks for and identifies any chain IDs
    that require concatenation within each MolID object.

    Parameters:
    -----------
    master_molID_class_list : list of MolID
        A list of MolID objects representing entity name and chain ID information from CIF(s).

    molIDConversion_list : list of MolIDConversion
        A list of MolIDConversion objects that provide the new chain ID mappings for each entity name.

    Returns:
    --------
    master_molID_class_list : list of MolID
        The updated list of MolID objects, with their chain ID mappings adjusted to the input
        MolIDConversion objects. Each MolID object is also checked for any chain ID concatenations.
    """
    molIDConvert_molID_chID_map = {}
    for molIDConversion in molIDConversion_list:
        molIDConvert_molID_chID_map[molIDConversion.molID] = molIDConversion.chID_list

    # Updates master_molID_class_list
    for molID_class in master_molID_class_list:
        for molID in molID_class.molID_chID:
            molID_class.add_chID_newchID_map(molID, molIDConvert_molID_chID_map[molID])

    # Determines which MolID objects now contain conflicts (concatenations)
    for molID_class in master_molID_class_list:
        molID_class.check_for_concatenations()
    return master_molID_class_list

# The class containing all the information about each file necessary to build
# a conversion template
class MolID(object):
    """
    A class to deal with and store entity name (MolID) information for a CIF,
    including mapping of chain IDs, concatenation orders, and completion statuses.

    Attributes:
    -----------
    file_name : str
        The name of the CIF.
    chID_newchID_map : dict
        A dictionary mapping original chain IDs (oldchID) to new chain IDs (newchID).
    molID_chID : dict
        A dictionary mapping entity names (molID) to their corresponding chain IDs (chID).
    concat_order : dict
        A dictionary tracking the order of concatenation for chain IDs, if applicable.
    complete_order : dict
        A dictionary indicating whether a chain ID has been fully processed (True) or not (False).

    Methods:
    --------
    add_chID_newchID_map(molID, newchID_list):
        Maps original chain IDs to new chain IDs for a given entity name.

    check_for_concatenations():
        Identifies chain IDs that need to be concatenated and updates the
        `concat_order` and `complete_order` dictionaries according to this.

    force_complete_order(chID, complete):
        Manually sets the completion status of a specific chain ID.

    update_concat_order(chID, neworder):
        Updates the concatenation order for a given chain ID, ensuring consistency with
        other chain IDs that share the same new chain ID.
    """
    def __init__(self, file_name, chID_newchID_map, molID_chID,
                 concat_order, complete_order):
        """
        Initializes the MolID class with the provided parameters.

        Parameters:
        -----------
        file_name : str
            The name of the CIF.
        chID_newchID_map : dict
            A dictionary to map original chain IDs (oldchID) to new chain IDs (newchID).
        molID_chID : dict
            A dictionary mapping entity names (molID) to their corresponding chain IDs (chID).
        concat_order : dict
            A dictionary to track the order of concatenation for chain IDs, if applicable.
        complete_order : dict
            A dictionary indicating whether a chain ID has been fully processed (True) or not (False).
        """
        self.file_name = file_name
        self.chID_newchID_map = chID_newchID_map
        self.molID_chID = molID_chID
        self.concat_order = concat_order
        self.complete_order = complete_order

    def add_chID_newchID_map(self, molID, newchID_list):
        """
        Maps original chain IDs to new chain IDs for a given entity name.

        Parameters:
        -----------
        molID : str
            The entity name for which the chain IDs are being mapped.
        newchID_list : list of str
            A list of new chain IDs corresponding to the original chain IDs.
        """
        N = 0
        for oldchID in self.molID_chID[molID]:
            self.chID_newchID_map[oldchID] = newchID_list[N]
            N += 1

    def check_for_concatenations(self):
        """
        Identifies chain IDs that need to be concatenated and updates the
        `concat_order` and `complete_order` dictionaries as so.
        """
        self.concat_order = {}
        # Figure out which newchIDs appear more than once which will be
        # concatenated
        usage = {}
        duplicates = {}
        for chID in self.chID_newchID_map:
            if self.chID_newchID_map[chID] not in usage:
                usage[self.chID_newchID_map[chID]] = 1
            else:
                duplicates[self.chID_newchID_map[chID]] = 1

        # Create a map of oldchID to newchID for those being concatenated
        usage = {}
        for chID in self.chID_newchID_map:
            if self.chID_newchID_map[chID] in duplicates:
                if self.chID_newchID_map[chID] in usage:
                    usage[self.chID_newchID_map[chID]] += 1
                else:
                    usage[self.chID_newchID_map[chID]] = 1
                self.concat_order[chID] = usage[self.chID_newchID_map[chID]]

        # Update complete_order
        for chID in self.chID_newchID_map:
            if chID in self.concat_order:
                self.complete_order[chID] = False
            else:
                self.complete_order[chID] = True


    # This will force the complete_order[chID] to the input "complete" which is either True or False
    def force_complete_order(self, chID, complete):
        """
        Sets the completion status of a specific chain ID.

        Parameters:
        -----------
        chID : str
            The original chain ID whose completion status is being set.
        complete : bool
            The completion status to be set for the chain ID (True or False).
        """
        self.complete_order[chID] = complete

    def update_concat_order(self, chID, neworder):
        """
        Updates the concatenation order for a given chain ID, ensuring consistency
        with other chain IDs that share the same new chain ID.

        Parameters:
        -----------
        chID : str
            The original chain ID whose concatenation order is being updated.
        neworder : int
            The new order value to assign to the chain ID.
        """
        otheroldchID = ""
        if chID in self.concat_order:
            newchID = self.chID_newchID_map[chID]
            for oldchID in self.chID_newchID_map:
                if (self.chID_newchID_map[oldchID] == newchID) and (oldchID != chID) and (self.concat_order[oldchID] == neworder):
                    otheroldchID = oldchID
            if (otheroldchID != ""):
                self.concat_order[otheroldchID] = self.concat_order[chID]
                self.concat_order[chID] = neworder

# END MolID class


# make_molID is the function that will grab all relevant information from each
# file input
def make_MolID(myfile): # NEVER CALLED
    """
    Process a file to extract entity names (MolID) and chain IDs (chID),
    and organize them into an instance of the `MolID` class.

    Function reads a file line by line, extracting entity name and
    chain ID information from lines containing "COMPND". The extracted
    data is then cleaned and structured into a dictionary mapping each
    entity name to its corresponding chain IDs. A new instance of the `MolID`
    class is created and returned.

    Parameters
    ----------
    myfile : file object
        The file object to be processed. It should contain lines with
        "COMPND" keywords where entity name and chID information are present.

    Returns
    -------
    my_molID_class : MolID
        An instance of the `MolID` class containing the file name, chain ID
        to new chain ID mapping, entity name (MolID) to chain ID mapping, concatenation
        order, and completion order.
    """
    file_name = myfile.name
    molID_list = []
    chID_list = []
    molID_chID = {}
    chID_newchID_map = {}
    concat_order = {}
    complete_order = {}

    for line in myfile:
        if "COMPND" in line:
            if "MOLECULE:" in line or "MOLID:" in line:
                line_list = line.rsplit(':')
                molID = line_list[1].strip()

                # Remove leading and trailing spaces
                molID = re.sub(r' \{0,\}$|^ \{0,\}', '', molID)

                # Remove any special characters
                molID = re.sub(':|;|’|“|”', '', molID)

                # Add single leading and trailing space and bookmark with ;
                molID = re.sub('^', '; ', molID)
                molID = re.sub('$', ' ;', molID)
                molID_list.append(molID.upper())

            if "CHAIN:" in line:
                line_list = line.rsplit(':')
                chID = line_list[1].strip()

                # Remove all spaces and ;
                chID = re.sub(' |;', '', chID)
                chID = chID.split(',')
                chID_list.append(chID)

    for i in range(len(chID_list)):
        if molID_list[i] not in molID_chID:
            molID_chID[molID_list[i]] = chID_list[i]
        else:
            for chid in chID_list[i]:
                    molID_chID[molID_list[i]].append(chid)

    my_molID_class = MolID(file_name, chID_newchID_map, molID_chID, concat_order, complete_order)
    my_molID_class.check_for_concatenations()
    return my_molID_class
# END make_MolID


# make_molID is the function that will grab all relevant information from each
# file input
def make_MolID_cif(myfile):
    """
    Extracts and maps entity names (MolID) to their chain IDs using the entity ID from a CIF,
    and returns an instance of the MolID class.

    Parameters:
    -----------
    myfile : file object
        A file object representing a CIF to be processed.

    Returns:
    --------
    my_molID_class : MolID
        An instance of the MolID class containing the extracted entity name's (MolIDs)
        with their chain IDs, and initialized but empty concatenation and completion order map.
    """
    file_name = myfile.name
    chID_newchID_map = {}
    molID_chID = {}
    concat_order = {}
    complete_order = {}

    mmcif_dict = MMCIF2Dict(myfile)

    # Link molID and chID using entity_id, creating a mapping between entity_id and auth_asym_id.
    # Store this mapping in entity_chIDlist_map

    entity_list = mmcif_dict['_atom_site.label_entity_id']
    chID_list = mmcif_dict['_atom_site.label_asym_id']
    entity_chIDlist_map = {}

    for i in range(len(entity_list)):
        if entity_list[i] not in entity_chIDlist_map:
            entity_chIDlist_map[entity_list[i]] = [chID_list[i]]
        else:
            if chID_list[i] not in entity_chIDlist_map[entity_list[i]]:
                entity_chIDlist_map[entity_list[i]].append(chID_list[i])

    entity_list = mmcif_dict['_entity.id']
    molID_list = mmcif_dict['_entity.pdbx_description']
    molID_list = [i.upper() for i in molID_list]

    for i in range(len(entity_list)):
        if entity_list[i] in entity_chIDlist_map:
            for chID in entity_chIDlist_map[entity_list[i]]:
                if molID_list[i] not in molID_chID:
                    molID_chID[molID_list[i]] = [chID]
                else:
                    molID_chID[molID_list[i]].append(chID)

    my_molID_class = MolID(file_name, chID_newchID_map, molID_chID, concat_order, complete_order)
    my_molID_class.check_for_concatenations()
    return my_molID_class
# END make_MolID_cif


# MolIDConversion Class. The idea here is to create an object that
# will be modified by the user in the first step to map a MolID (entity name) to a given
# chain ID with the correct number of occurrences. This is the master list of
# molID's (entity names)
class MolIDConversion(object):
    """
    A class that handles the conversion of entity names (MolIDs) to chain IDs,
    tracking the completeness of the conversion based on the expected number
    of occurrences.

    Attributes:
    -----------
    molID : str
        The entity name (MolID) being mapped.
    chID_list : list of str
        A list of chain IDs associated with the entity name.
    occur : int
        The expected number of chain ID occurrences for the entity name (MolID).
    complete : bool
        A boolean indicating whether the conversion is complete.

    Methods:
    --------
    check_for_completeness():
        Checks if the number of chain IDs in the list meets or exceeds the expected
        occurrences and updates the completeness status.

    add_chID(chID):
        Adds a single chain ID to the `chID_list`.

    add_chID_list(chID_list):
        Adds multiple chain IDs to the `chID_list`.

    remove_chID_list(chID_list):
        Removes specified chain IDs from the `chID_list` if they exist.
    """
    molID = ""
    chID_list = []
    occur = 0
    complete = bool

    def __init__(self, molID, chID_list, occur, complete):
        """
        Initializes the MolIDConversion class with the provided parameters.

        Parameters:
        -----------
        molID : str
            The entity name (MolID) being mapped.
        chID_list : list of str
            A list of chain IDs associated with the entity name.
        occur : int
            The number of chain ID occurrences for the entity name.
        complete : bool
            A boolean indicating whether the conversion is 'complete'.
        """
        self.molID = molID
        self.chID_list = chID_list
        self.occur = occur
        self.complete = complete

    def check_for_completeness(self):
        """
        Checks if the number of chain IDs in the list matches or goes over the expected
        occurrences and updates the completeness status.

        If the number of chain IDs in `chID_list` is greater than or equal to `occur` (occurrences of the entity
        in the cif), `complete` (all chain IDs have been assigned with a chain ID assignment) is set to True. Otherwise
        , it is set to False.
        """
        if (len(self.chID_list) >= self.occur):
            self.complete = True
        else:
            self.complete = False

    def add_chID(self, chID): # NEVER CALLED
        """
        Adds a chain ID to the `chID_list`.

        Parameters:
        -----------
        chID : str
            The chain ID to be added to the list.
        """
        self.chID_list.append(chID)

    def add_chID_list(self, chID_list):
        """
        Adds multiple chain IDs to the `chID_list`.

        Parameters:
        -----------
        chID_list : list of str
            A list of chain IDs to be added to the list.
        """
        for chID in chID_list:
            self.chID_list.append(chID)

    def remove_chID_list(self, chID_list):
        """
        Removes specified chain IDs from the `chID_list` if they exist.

        Parameters:
        -----------
        chID_list : list of str
            A list of chain IDs to be removed from the list.
        """
        for chID in chID_list:
            if chID in self.chID_list:
                self.chID_list.remove(chID)


def initial_MolIDConversion(molID, occur):
    """
    Initializes a MolIDConversion object with the entity name (MolID) and expected number of chain ID occurrences.

    Parameters:
    -----------
    molID : str
        The entity name being mapped.
    occur : int
        The expected number of chain ID occurrences for the entity name.

    Returns:
    --------
    my_molID_conversion : MolIDConversion
        An initialized MolIDConversion object with the provided entity name and expected occurrences.
    """
    chID_list = []
    complete = False
    my_molID_conversion = MolIDConversion(molID, chID_list, occur, complete)
    return my_molID_conversion


def add_chID_list(MolIDConversion, chID_list):
    """
    Adds a list of chain IDs to a MolIDConversion object, making sure that there's no duplicates.

    Parameters:
    -----------
    MolIDConversion : MolIDConversion
        An instance of the MolIDConversion class to which the chain IDs will be added.

    chID_list : list of str
        A list of chain IDs to be added to the MolIDConversion object's `chID_list`.

    Returns:
    --------
    MolIDConversion : MolIDConversion
        The updated MolIDConversion object with the new chain IDs added.
    """
    for chID in chID_list:
        if chID not in MolIDConversion.chID_list:
            MolIDConversion.chID_list.append(chID)
    return MolIDConversion
#

# Compile information from each file to create a master map of molID to max
# number of occurrences
def CreateMasterUniqueMolIDMap(molID_class_list):
    """
    Compiles information from CIF(s) to create a master map of entity names to the max number of occurrences
    of their associated chain ID's.

    Parameters:
    -----------
    molID_class_list : list
        A list of `MolID` class instances, each representing data extracted from a CIF
        (mapping of entity name to chain ID).

    Returns:
    --------
    unique_molID_map : dict
        A dictionary where keys are `molID`(entity name) strings, and values are the maximum number
        of occurrences of their associated chain IDs across all CIF(s).
    """
    #essentially counts the amount of chain id's each entity name has
    #it iterates over every cif file
    unique_molID_map = {}
    for my_molID_class in molID_class_list:
        for molID in my_molID_class.molID_chID:
            if (unique_molID_map.get(molID) is None):
                unique_molID_map[molID] = \
                    len(my_molID_class.molID_chID[molID])
            elif (unique_molID_map.get(molID) <
                    len(my_molID_class.molID_chID[molID])):
                unique_molID_map[molID] = len(my_molID_class.molID_chID[molID])
    return unique_molID_map

# FAPA APRIL 2024 V2 STARTS

# Compile information from each file to create a master map of molID to
# number of occurrences
def CreateMasterUniqueMolIDOccursLIST(molID_class_list):
    """
    Creates a dictionary that maps each entity name (MolID) to a list of occurrence counts.

    Parameters:
    -----------
        molID_class_list : list
        A list of objects, each containing a `molID_chID` dictionary that
        maps entity names to chain IDs.

    Returns:
    --------
        unique_molID_map_list : dict
        A dictionary where the keys are entity names (MolIDs) and the values are lists containing the number of
        chain IDs associated with each entity name in different CIF(s).
    """
    unique_molID_map_list = {}
    for my_molID_class in molID_class_list:
        for molID in my_molID_class.molID_chID:
            if (unique_molID_map_list.get(molID) is None):
                unique_molID_map_list[molID] = \
                    [len(my_molID_class.molID_chID[molID])]
            else:
                unique_molID_map_list[molID].append(len(my_molID_class.molID_chID[molID]))
    return unique_molID_map_list

def CreateMasterUniqueMolIDinitialChainIDsLIST(molID_class_list):
    """
    Creates a dictionary that maps each entity name (MolID) to a list of initial chain IDs.

    Parameters
    ----------
    molID_class_list : list
        A list of objects, each containing a `molID_chID` dictionary that maps MolIDs to chain IDs.

    Returns
    -------
    unique_molID_ChainIDs_map_list : dict
        A dictionary where the keys are entity names  and the values are lists of lists,
        each containing the chain IDs associated with that entity name in different CIF(s).
    """
    unique_molID_ChainIDs_map_list = {}
    for my_molID_class in molID_class_list:
        for molID in my_molID_class.molID_chID:
            if (unique_molID_ChainIDs_map_list.get(molID) is None):
                unique_molID_ChainIDs_map_list[molID] = \
                    [my_molID_class.molID_chID[molID]]
            else:
                unique_molID_ChainIDs_map_list[molID].append([my_molID_class.molID_chID[molID]])
    return unique_molID_ChainIDs_map_list


# FAPA APRIL 2024 V2 END

# FAPA MARCH 2024
# This is a test to create a "track-changes" feature
# We will keep track of the structures that contain
# each entity :D let's see!
def CreateMasterUniqueMolIDMapWithFileName(molID_class_list):
    """
    Creates a dictionary that maps each unique entity name (MolID) to a list
    of CIF(s) containing that entity name.

    Parameters:
    -----------
        molID_class_list : list
            A list of objects, each containing a `molID_chID` dictionary that maps
            entity names to chain IDs, and a `file_name` attribute representing the name of the CIF.

    Returns:
    --------
        files_contain_molID : dict
        A dictionary where the keys are entity names and the values are lists of CIF names that
         contain the entity name.
    """
    files_contain_molID = {} #this will be a dictionary of entity:[list of files]
    for my_molID_class in molID_class_list:
        for molID in my_molID_class.molID_chID:
            if molID not in files_contain_molID:
                files_contain_molID[molID] = [my_molID_class.file_name]
            else:
                files_contain_molID[molID].append(my_molID_class.file_name)
    return files_contain_molID

def show_full_conversion_and_file_list_by_number_chains(MolID_ChainID_dict_of_lists,current_MolID_class_list,current_MolID_file_list,MolID_occur_dict_of_lists,target_dir,write_csv=True):
    """
    Display and optionally save a mapping of old chain IDs to new chain IDs, including
    associated file lists and the number of occurrences of each entity name.

    Parameters
    ----------
    MolID_ChainID_dict_of_lists : dict
        A dictionary where each key is a entity name, and the corresponding value is a list
        of old chain IDs for that entity name.
    current_MolID_class_list : list
        A list of `MolIDConversion` class instances, each representing the mapping between old
        and new chain IDs for a entity name.
    current_MolID_file_list : dict
        A dictionary where each key is a entity name, and the corresponding value is a list
        of file paths associated with that entity name.
    MolID_occur_dict_of_lists : dict
        A dictionary where each key is a entity name, and the corresponding value is a list
        of integers representing the number of occurrences (chains) in the corresponding file.
    target_dir : str
        The directory path where the CSV file will be saved if `write_csv` is `True`.
    write_csv : bool, optional
        If `True`, the function will save the output to a CSV file. Default is `True`.

    Returns
    -------
    None
    """
    if write_csv:
        with open(f'{target_dir}/OldChainID_NewChainID_numbered_MolID_Files_Map.csv', 'w') as fout:
            fout.write('OldChainID(label_asym_id):NewChainID:Entity:NumberOfFiles:Files\n')

    for molIDConversion in current_MolID_class_list:
        filelist = [x.split("/")[-1] for x in current_MolID_file_list[molIDConversion.molID]]
        counter=0
        for file_name in filelist:
            nchains=MolID_occur_dict_of_lists[molIDConversion.molID][counter]
            if nchains > 5:
                molIDCon_chID_list_forPrint = str(nchains)+"x"+str(molIDConversion.chID_list[0])
            else:
                molIDCon_chID_list_forPrint = re.sub(r"\[|\]| |'", '', str(molIDConversion.chID_list[:nchains]))

            chainds_for_print = re.sub(r"\[|\]| |'", '',str(MolID_ChainID_dict_of_lists[molIDConversion.molID][counter]))

            print(chainds_for_print + ":" + molIDCon_chID_list_forPrint+":"+molIDConversion.molID+":"+str(len(filelist))+":"+file_name)
            if write_csv:
                with open(f'{target_dir}/OldChainID_NewChainID_numbered_MolID_Files_Map.csv', 'a') as fout:
                    fout.write(f'{chainds_for_print}:{molIDCon_chID_list_forPrint}:{molIDConversion.molID}:{str(len(filelist))}:{file_name}\n')
            counter+=1

    print("\n")
#
#
# Read input file function
def read_input_file(input_cnv_file):
    """
    Reads a user input file to create a mapping of entity names (MolIDs) to their chain IDs.

    Parameters:
    -----------
    input_cnv_file : str
        The path to the input file containing conversion data. Each line should have
        an entity name (MolID) followed by chain IDs, separated by colons and commas.

    Returns:
    --------
    user_molID_chID_map : dict
        A dictionary where keys are entity names (MolIDs) and values are lists of chain IDs
        associated with each entity name.
    """
    user_molID_chID_map = {}
    user_molID_list_v = []
    user_chID_list_v = []

    # Read user input file and create user_molID_chID_map
    if (input_cnv_file != ""):
        if (os.path.isfile(input_cnv_file) is True):
            # Reading of input file
            with open(input_cnv_file) as mycnv:
                for line in mycnv:
                    line = line.strip()
                    my_line = line.split(':')
                    molID = my_line[0]
                    chID = my_line[1].split(',')
                    user_molID_list_v.append(molID)
                    user_chID_list_v.append(chID)

            #links the entity name to the chain ID
            for molID_i, chID_i in zip(user_molID_list_v, user_chID_list_v):
                if molID_i not in user_molID_chID_map:
                    user_molID_chID_map[molID_i] = chID_i
                else:
                    for this_chID in chID_i:
                        if this_chID not in user_molID_chID_map[molID_i]:
                            user_molID_chID_map[molID_i].append(this_chID)

        return user_molID_chID_map
# End Read user input file function


####################################
# INTERACTIVE CONVERSION FUNCTIONS #
####################################

def show_full_conversion(current_list, step='conversion'):
    """
    Displays the full conversion or concatenation information of entity names and their chain IDs.

    Parameters:
    -----------
    current_list : list
        - A list of `molIDConversion` objects if `step` is 'conversion'.
        - A list of `molID` class objects if `step` is 'concatenation'.
    step : str, optional
        Specifies the mode of operation. Default is 'conversion'.

    Returns:
    --------
    None
    """
    if(step=='conversion'):
        for molIDConversion in current_list:
            molIDCon_chID_list_forPrint = re.sub(r"\[|\]| |'", '', str(molIDConversion.chID_list))
            print(molIDConversion.molID+":"+molIDCon_chID_list_forPrint)
    elif(step=='concatenation'):
        for molID_class in current_list:
            for molID in molID_class.molID_chID:
                for chID in molID_class.molID_chID[molID]:
                    if chID in molID_class.concat_order:
                        print(molID_class.file_name + ":" + molID + ":"
                              + chID + ":" +
                              molID_class.chID_newchID_map[chID] + ":"
                              + str(molID_class.concat_order[chID]))
                    else:
                        print(molID_class.file_name + ":" + molID + ":"
                              + chID + ":" +
                              molID_class.chID_newchID_map[chID] + ":" +
                              str(0))

def show_unassigned_conversion(current_list, step='conversion'):
    """
    Displays unassigned entity name conversions or concatenations.

    Parameters:
    -----------
    current_list : list
        A list of objects, either `MolIDConversion` or `MolID_class`, depending on the step.
        Each object in the list will be checked for completeness based on the specified step.

    step : str, optional
        Specifies the type of check to perform. It can be either 'conversion' or 'concatenation'.
        - 'conversion': Checks for incomplete entity name conversions.
        - 'concatenation': Checks for incomplete entity name concatenations.
        """
    if(step=='conversion'):
        counter1=0
        counter2=0
        for molIDConversion in current_list:
            molIDConversion.check_for_completeness()
            if molIDConversion.complete is False:
                counter1+=1
                counter2=counter2+molIDConversion.occur
                molIDCon_chID_list_forPrint = re.sub(r"\[|\]| |'", '', str(molIDConversion.chID_list))
                #
                print(str(molIDConversion.occur)+":"+str(molIDConversion.molID)+":"+molIDCon_chID_list_forPrint)
        print("You need to accept %s entity conversions" % counter1)
        print("You need to accept %s total chain conversions" % counter2)

    elif(step=='concatenation'):
        counter=0
        for molID_class in current_list:
            for molID in molID_class.molID_chID:
                for chID in molID_class.molID_chID[molID]:
                    if molID_class.complete_order[chID] is False:
                        counter+=1
                        print(molID_class.file_name + ":" + molID + ":"
                              + chID + ":" +
                              molID_class.chID_newchID_map[chID] + ":"
                              + str(molID_class.concat_order[chID]))
        print("You need to accept %s concatenations" % counter )


def return_unassigned_conversion(current_list, step='conversion'):
    """
    Identifies and prints or stores unassigned conversions or concatenations.

    Function processes a list of conversion or concatenation objects to
    determine if any chain ID assignments are incomplete. It either prints
    the incomplete conversions or stores the unassigned concatenations in a list.

    Parameters:
    -----------
    current_list : list
        A list of either `molIDConversion` objects (for conversion) or
        `molID_class` objects (for concatenation), each containing information
        about entity name (MolID) and associated chain IDs.

    step : str, optional
        Specifies the operation mode. The default is 'conversion'.
        - 'conversion' : Checks and prints incomplete conversions.
        - 'concatenation' : Checks and stores incomplete concatenations in a list.

    Returns:
    --------
    unassigned : list
        A list of unassigned concatenations if `step` is 'concatenation'. Each
        entry in the list is a string containing the file name, MolID, chain ID,
        new chain ID, and the concatenation order. If `step` is 'conversion',
        this list will be empty.
    """
    unassigned = []

    if(step=='conversion'):
        for molIDConversion in current_list:
            molIDConversion.check_for_completeness()
            if molIDConversion.complete is False:
                molIDCon_chID_list_forPrint = re.sub(r"\[|\]| |'", '', str(molIDConversion.chID))#FAPA
                print(str(molIDConversion.occur)+":"+str(molIDConversion.molID)+":"+molIDCon_chID_list_forPrint)
    elif(step=='concatenation'):
        for molID_class in current_list:
            for molID in molID_class.molID_chID:
                for chID in molID_class.molID_chID[molID]:
                    if molID_class.complete_order[chID] is False:
                        unassigned.append(molID_class.file_name + ":" + molID + ":"
                              + chID + ":" +
                              molID_class.chID_newchID_map[chID] + ":"
                              + str(molID_class.concat_order[chID]))
    return unassigned

def add_user_conversion(molIDConversion_list):
    """
    Updates a list of MolIDConversion objects with chain IDs from a user-provided input file
    and checks each for completeness.

    This function takes a list of MolIDConversion objects, reads the user input file
    that has the entity names and their chain IDs, then adds the chain IDs to the corresponding
    MolIDConversion objects. After updating, the function checks each MolIDConversion object
    for completeness based on the updated chain ID lists.

    Parameters:
    -----------
    molIDConversion_list : list
        molIDConversion objects where each object contains an entity name and their chain IDs

    Returns:
    -----------
    molIDConversion_list : list
        Updated list of molIDConversion objects after adding the chain IDs from the users input file
    """
    input_cnv_file = input('Conversion File: ')
    user_molID_chID_map =  read_input_file(input_cnv_file)
    for molIDConversion in molIDConversion_list:
        # This is currently strict inclusion but perhaps should be except out
        # of convenience to the user
        for key in user_molID_chID_map:
            if (str(key)==str(molIDConversion.molID)):
                molIDConversion.add_chID_list(user_molID_chID_map[key])
        molIDConversion.check_for_completeness()
    # !! molIDConversion_list has been updated
    return molIDConversion_list

def edit_conversion_interface(molIDConversion_list, action='add'): #FAPA: IT WAS A TYPO!!! -_-
    """
    Interface for interacting with and modifying molIDConversion objects in the list.

    This function allows users to search for molID objects, narrow down the search results, and modify chain IDs based
    on the action ('add' or 'remove').

    Parameters:
    -----------
    molIDConversion_list : list
        A list of molIDConversion objects containing entity names (molID) and their associated chain IDs.

    action : str, optional
        The action to perform on the chain IDs. It determines whether chain IDs should be added or removed.
        Default is 'add'.

    Returns:
    --------
    molIDConversion_list : list
        The updated list of molIDConversion objects after performing the specified action.
    """
    search_term = input('MolID search term: ')
    molIDConversion_list, search_molIDConversion_list = search_conversion(molIDConversion_list, search_term) #FAPA removed molidutils.
    input_submenu = 0
    while (input_submenu != "DONE"):
        print("    1) Further narrow down search results")
        if(action=='add'):
            print("    2) Add chain ID to conversion templates")
        elif(action=='remove'):
            print("    2) Remove chain ID from conversion templates")
        input_submenu = input('Option Number: ')
        if (input_submenu == "QUIT"):
            input_submenu = "DONE"
        if (input_submenu == "1"):
            search_term = input('MolID search term: ')
            molIDConversion_list, search_molIDConversion_list = search_again_conversion(molIDConversion_list, search_molIDConversion_list, search_term)
        elif (input_submenu == "2"):
            if(action=='add'):
                print("""Enter new chain IDs, comma separated, no spaces""")
            elif(action=='remove'):
                print("""Enter chain ID to be removed, comma separated, no spaces""")
            chID_list = input('Chain IDs: ')
            if (chID_list == "") or (chID_list == "QUIT"):
                pass
            else:
                molIDConversion_list, search_molIDConversion_list = edit_chain_conversion(molIDConversion_list, search_molIDConversion_list, chID_list, action=action)
            input_submenu = "DONE"
    return molIDConversion_list

def search_conversion(molIDConversion_list, search_term):
    """
    Search for an entity name in the molIDConversion list, print its details, and remove it from the list.

    This function checks whether the entity name (search term) exists in the list of
    molIDConversion objects. If a match is found, the corresponding molIDConversion object is appended
    to a new list, its details (entity name and chain IDs) are printed, and it is removed from the original
    molIDConversion_list.

    Parameters:
    -----------
    molIDConversion_list : list
        A list of molIDConversion objects consisting of the entity name with their chain IDs
    search_term : str
        The entity name to search for in the molIDConversion list.

    Returns:
    --------
    molIDConversion_list : list
        The updated list with the searched molIDConversion objects removed.
    search_molIDConversion_list : list
        A list of the molIDConversion objects that matched the search term.
    """
    search_molIDConversion_list = []
    # checks to see that the entity name they looked up is in the list of entity names we saved,
    #
    for molIDConversion in molIDConversion_list:
        if search_term in molIDConversion.molID:
            search_molIDConversion_list.append(molIDConversion)
            molIDCon_chID_list_forPrint = re.sub(r"\[|\]| |'", '', str(molIDConversion.chID_list))
            print(molIDConversion.molID+":"+molIDCon_chID_list_forPrint)
    for molIDConversion in search_molIDConversion_list:
        molIDConversion_list.remove(molIDConversion)
    return molIDConversion_list, search_molIDConversion_list

def search_again_conversion(molIDConversion_list, search_molIDConversion_list, search_term):
    """
    Searches for an entity name (molID) in the search_molIDConversion_list and updates both lists based on the search
    results.

    If search term is found in search_molIDConversion_list, the corresponding molIDConversion object is added to a
    new list.
    If not found, the molIDConversion object is added to molIDConversion_list.

    Parameters:
    -----------
    molIDConversion_list : list
        A list of molIDConversion objects, each with an entity name (molID) with its chain IDs.

    search_molIDConversion_list : list
        A list of molIDConversion objects that have been previously searched or filtered out from molIDConversion_list.

    search_term : str
        The entity name (molID) that the user wants to search for within the search_molIDConversion_list.

    Returns:
    --------
    molIDConversion_list : list
        The updated list of molIDConversion objects, with the search term entity either still removed or added back,
        depending on the search result.

    search_molIDConversion_list : list
        The updated list containing only the molIDConversion objects that match the search term.
    """
    new_search_molIDConversion_list = []
    for molIDConversion in search_molIDConversion_list:
        if search_term in molIDConversion.molID:
            new_search_molIDConversion_list.append(molIDConversion)
            molIDCon_chID_list_forPrint = re.sub(r"\[|\]| |'", '', str(molIDConversion.chID_list))
            print(molIDConversion.molID+":"+molIDCon_chID_list_forPrint)

        else:
            molIDConversion_list.append(molIDConversion)

    search_molIDConversion_list = new_search_molIDConversion_list
    return molIDConversion_list, search_molIDConversion_list

def edit_chain_conversion(molIDConversion_list,search_molIDConversion_list, chID_list, action='add'):
    """
    Edits the chain ID list of molIDConversion objects based on the action and updates the main list.

    This function performs an addition or removal of chain IDs to/from the `search_molIDConversion_list` based on the
    provided action. After modifying the chain IDs, the function adds the updated `molIDConversion` objects back to the
    `molIDConversion_list` and checks their completeness.

    Parameters:
    -----------
    molIDConversion_list : list
        A list of molIDConversion objects that represents the complete list of entity names and their associated chain
        IDs.

    search_molIDConversion_list : list
        A list of molIDConversion objects that are currently being edited. These are the entities for which chain IDs
        will be modified.

    chID_list : str
        A comma-separated string of chain IDs to be added or removed from the `search_molIDConversion_list` objects.

    action : str, optional
        The action to perform on the chain IDs. Options are 'add' to add chain IDs and 'remove' to remove chain IDs.
        Default is 'add'.

    Returns:
    --------
    molIDConversion_list : list
        The updated list of molIDConversion objects, including the modified objects from `search_molIDConversion_list`.

    search_molIDConversion_list : list
        The updated list of molIDConversion objects with modified chain IDs, reflecting the changes made by the action.
    """
    chID_list = chID_list.split(',')
    # Performing chain addition to searched terms
    for molIDConversion in search_molIDConversion_list:
        if(action=='add'):
            molIDConversion.add_chID_list(chID_list)
        elif(action=='remove'):
            molIDConversion.remove_chID_list(chID_list)
    # Adding modified molIDConversions back to molIDConversion_list
    for molIDConversion in search_molIDConversion_list:
        molIDConversion.check_for_completeness()
        molIDConversion_list.append(molIDConversion)
    return molIDConversion_list, search_molIDConversion_list

def edit_conversion_manual(molIDConversion_list):
    """
    Manually edit chain IDs for each molIDConversion object in the list.

    This function asks user to enter chain IDs for each entity name (molIDConversion object).
    Each entity name is updated with the new chain IDs and checks their completeness.

    Parameters:
    -----------
    molIDConversion_list : list
        A list of molIDConversion objects, each representing an entity name (molID) with its associated chain IDs.

    Returns:
    --------
    molIDConversion_list : list
        The updated list of molIDConversion objects, with each object modified based on the user provided chain IDs.
    """
    print("Enter chain IDs for each of the following MolID.")
    print("Comma separated, no spaces")
    for molIDConversion in molIDConversion_list:
        chID_list = input(molIDConversion.molID+":")
        if (chID_list == "") or (chID_list == "QUIT"):
            pass
        else:
            chID_list = chID_list.split(',')
        molIDConversion.add_chID_list(chID_list)
        molIDConversion.check_for_completeness()
    return molIDConversion_list

def check_complete(molIDConversion_list):
    """
    Counts the number of entity names with incomplete chain ID assignments and indicates if all are complete.

    This function iterates through a list of molIDConversion objects and counts how many of them have incomplete chain
    ID assignments.

    Parameters:
    -----------
    molIDConversion_list : list
        A list of molIDConversion objects, each representing an entity name with chain ID assignments.

    Returns:
    --------
    input_menu_complete : str
        A string indicating the completeness of chain ID assignments:
        - "0" if there are any incomplete assignments.
        - "1" if all assignments are complete.
    """
    num_unassigned = 0
    for molIDConversion in molIDConversion_list:
        if molIDConversion.complete is False:
            num_unassigned += 1
    if (num_unassigned > 0):
        input_menu_complete = "0"
    else:
        input_menu_complete = "1"
    return input_menu_complete

#######################################
# INTERACTIVE CONCATENATION FUNCTIONS #
#######################################

def problem_counter(master_molID_class_list):
    """
    Counts the number of incomplete entries in the master list of molID classes.

    This function iterates through each molID class in the provided list and counts how many chain IDs in the
    `complete_order` attribute are marked as `False`, indicating incomplete entries

    Parameters:
    -----------
    master_molID_class_list : list
        A list of molID class objects, where each object contains a `complete_order` dictionary that maps chain IDs to
        completion statuses.

    Returns:
    --------
    count_problems : int
        The number of incomplete chain ID entries (marked as `False`) across all molID class objects in the list.
    """
    count_problems = 0
    for molID_class in master_molID_class_list:
        for chID in molID_class.complete_order:
            if molID_class.complete_order[chID] is False:
                count_problems += 1
    return count_problems

def edit_concatenation_interface(master_molID_class_list, new_order=None, action='accept'):
    """
    Manages interface for editing or accepting chain ID concatenation.

    Function provides an interactive interface for users to search for
    entity names, update chain IDs, modify concatenation order, or accept
    proposed concatenations. Depending on the 'action', the user can
    try new concatenation scenarios, update the existing concatenation order,
    or accept the proposed plan.

    Parameters:
    -----------
    master_molID_class_list : list
        A list of `molID_class` objects representing different entity names
        and their associated chain IDs.

    new_order : str, optional
        The new chain ID or concatenation order provided by the user, depending
        on the operation mode. Defaults is None.

    action : str, optional
        Specifies the type of operation:
        - 'try' : Allows user to test new chain IDs.
        - 'update' : Allows user to update the concatenation order.
        - 'accept' : Allows user to accept the current concatenation plan.
        Default 'action' is 'accept'.

    Returns:
    --------
    master_molID_class_list : list
        The updated list of `molID_class` objects after applying the user's
        modifications.

    new_order : str
        The last chain ID or concatenation order input by the user. If no new
        input was provided, it returns the initial value of `new_order`.
    """
    concat_submenu = 0
    while(concat_submenu != "QUIT"):
        search_term, concat_submenu = get_search_term(concat_submenu)
        if concat_submenu != "QUIT":
            found_molID_class_chID_map, molID_class_been_copied = search_chains(master_molID_class_list, search_term)
        while (concat_submenu != "QUIT"):
            if(action=='try'):
                print("""Select one of the following options to proceed:
                         1) Perform new search
                         2) Update new chain ID
                      """)
            elif(action=='update'):
                print("""Select one of the following options to proceed:
                         1) Perform new search
                         2) Update concatenation order
                      """)
            elif(action=='accept'):
                print("""Select one of the following options to proceed:
                         1) Perform new search
                         2) Accept planned concatenation
                      """)
            concat_submenu = input('Option Number: ')
            if (concat_submenu == "1"):
                break
            elif (concat_submenu == "2"):
                if(action=='try'):
                    new_order = input('New Chain ID: ')
                elif(action=='update'):
                    new_order = input('New concatenation order: ')
                found_molID_class_chID_map = edit_chain_order(found_molID_class_chID_map, new_order, action=action)
                if(action=='try'):
                    print_conflicts(found_molID_class_chID_map)
                    print("Would you like to accept or deny these changes?")
                    while (concat_submenu != "QUIT"):
                        concat_submenu = input('Enter ACCEPT or DENY: ')
                        if (concat_submenu == "ACCEPT"):
                            master_molID_class_list = accept_newchain(master_molID_class_list, found_molID_class_chID_map)
                        concat_submenu = "QUIT"

                else:
                    master_molID_class_list = accept_newchain(master_molID_class_list, found_molID_class_chID_map)
                    concat_submenu = "QUIT"
    return master_molID_class_list, new_order


def list_accept_concatenations(master_molID_class_list, search_term, new_order=None, action='accept'):
    """
    Manages and resolves chain concatenations in MolID classes within a master list.

    Allows for user to search for chains in the MolID classes, review planned concatenations,
    and either accept, deny, or modify these concatenations based on user input.

    Parameters:
    -----------
    master_molID_class_list : list of MolID
        A list containing instances of MolID classes representing entity names and their associated chains.
    search_term : str
        The term used to search for specific chains within the MolID classes. This is in the format
        "MolID:ChainID".
    new_order : str, optional
        The new chain ID or concatenation order, depending on the action, that the user wants to apply. Default is None.
    action : str, optional
        The action to perform on the chains. It can be:
        - 'try': Allows user to try out a new chain ID.
        - 'update': Allows user to update the concatenation order.
        - 'accept': Allows user to accept the planned concatenation. Default is 'accept'.

    Returns:
    --------
    master_molID_class_list : list of MolID
        The updated master list of MolID classes after applying the selected changes.
    new_order : str
        The new chain ID or concatenation order applied, if any.
    """
    concat_submenu = 0
    while(concat_submenu != "QUIT"):
        search_term = search_term.split(":")

        if concat_submenu != "QUIT":
            found_molID_class_chID_map, molID_class_been_copied = search_chains(master_molID_class_list, search_term)
        while (concat_submenu != "QUIT"):
            if(action=='try'):
                print("""Select one of the following options to proceed:
                         1) Perform new search
                         2) Update new chain ID
                      """)
            elif(action=='update'):
                print("""Select one of the following options to proceed:
                         1) Perform new search
                         2) Update concatenation order
                      """)
            elif(action=='accept'):
                print("""Select one of the following options to proceed:
                         1) Perform new search
                         2) Accept planned concatenation
                      """)
            concat_submenu = input('Option Number: ')
            if (concat_submenu == "1"):
                break
            elif (concat_submenu == "2"):
                if(action=='try'):
                    new_order = input('New Chain ID: ')
                elif(action=='update'):
                    new_order = input('New concatenation order: ')
                found_molID_class_chID_map = edit_chain_order(found_molID_class_chID_map, new_order, action=action)
                if(action=='try'):
                    print_conflicts(found_molID_class_chID_map)
                    print("Would you like to accept or deny these changes?")
                    while (concat_submenu != "QUIT"):
                        concat_submenu = input('Enter ACCEPT or DENY: ')
                        if (concat_submenu == "ACCEPT"):
                            master_molID_class_list = accept_newchain(master_molID_class_list, found_molID_class_chID_map)
                        concat_submenu = "QUIT"

                else:

                    master_molID_class_list = accept_newchain(master_molID_class_list, found_molID_class_chID_map)

                    concat_submenu = "QUIT"
    return master_molID_class_list, new_order


def list_accept_concatenations_auto(master_molID_class_list, search_term, new_order=None, action='accept'):
    """
    Automate the process of handling concatenation options for entity names based on user input.

    Function allows uses to perform a new search, update chain IDs, or accept planned concatenations
    based on the search term. It asks user to accept or update concatenation details and
    modifies the entity name class list as so.

    Parameters
    ----------
    master_molID_class_list : list
        A list of `MolID` class instances, each containing entity name and chain ID information.
    search_term : str
        A search term formatted as 'File:MolID:OldChain:NewChain:ConcatOrder' used to filter the entity names
        and chain IDs for concatenation operations.
    new_order : str, optional
        A string representing a new concatenation order or chain ID, depending on the `action`. Default is None.
    action : str, optional
        The action to be performed, which can be one of 'try', 'update', or 'accept'.
        - 'try': Attempts to update the chain ID or concatenation order.
        - 'update': Updates the concatenation order.
        - 'accept': Accepts the planned concatenation changes.
        Default is 'accept'.

    Returns
    -------
    master_molID_class_list : list
        The updated list of `MolID` class instances after processing the concatenation options.
    new_order : str
        The updated new order or chain ID provided by the user, if applicable.
    """
    concat_submenu = 0
    while(concat_submenu != "QUIT"):
        search_term = search_term.split(":")

        if concat_submenu != "QUIT":
            found_molID_class_chID_map, molID_class_been_copied = search_chains(master_molID_class_list, search_term)
        while (concat_submenu != "QUIT"):
            if(action=='try'):
                print("""Select one of the following options to proceed:
                         1) Perform new search
                         2) Update new chain ID
                      """)
            elif(action=='update'):
                print("""Select one of the following options to proceed:
                         1) Perform new search
                         2) Update concatenation order
                      """)
            elif(action=='accept'):
                print("""Select one of the following options to proceed:
                         1) Perform new search
                         2) Accept planned concatenation
                      """)

            concat_submenu = "2"
            if (concat_submenu == "1"):
                break
            elif (concat_submenu == "2"):
                if(action=='try'):
                    new_order = input('New Chain ID: ')
                elif(action=='update'):
                    new_order = input('New concatenation order: ')
                found_molID_class_chID_map = edit_chain_order(found_molID_class_chID_map, new_order, action=action)
                if(action=='try'):
                    print_conflicts(found_molID_class_chID_map)
                    print("Would you like to accept or deny these changes?")
                    while (concat_submenu != "QUIT"):
                        concat_submenu = input('Enter ACCEPT or DENY: ')
                        if (concat_submenu == "ACCEPT"):
                            master_molID_class_list = accept_newchain(master_molID_class_list, found_molID_class_chID_map)
                        concat_submenu = "QUIT"

                else:

                    master_molID_class_list = accept_newchain(master_molID_class_list, found_molID_class_chID_map)

                    concat_submenu = "QUIT"
    return master_molID_class_list, new_order


def get_search_term(value):
    """
    Prompts the user to input a search term in a specific format until the input is valid or the user chooses to quit.

    Parameters:
    -----------
    value : str
        A string representing an initial value passed to the function, which will be updated based on user input.

    Returns:
    --------
    search_term : list of str
        A list containing the components of the search term split by colons.
    value : str
        Updated value based on user input. If input is "QUIT", the value is set to "QUIT".
    """
    search_ok = ""
    while (search_ok != "OK"):
        print("Search format - File:MolID:OldChain:NewChain:ConcatOrder")
        search_term = input('Search ')
        if search_term == "QUIT":
            search_ok = "OK"
            value = "QUIT"
        search_term = search_term.split(":")
        if (len(search_term) > 4):
            search_ok = "OK"
    return search_term, value


def get_concat_line_term(value): # NEVER CALLED
    """
    Prompt the user to input a concatenation search term in the format
    'File:MolID:OldChain:NewChain:ConcatOrder' and return the processed search term.

    Parameters
    ----------
    value : str
        An initial value that may be updated to "QUIT" if the user decides to exit the loop.

    Returns
    -------
    search_term : list
        A list of strings representing the search term components:
        [File, MolID, OldChain, NewChain, ConcatOrder].
    value : str
        The updated value, which will be "QUIT" if the user exits the loop.
    """
    search_ok = ""
    while (search_ok != "OK"):
        print("concat format - File:MolID:OldChain:NewChain:ConcatOrder")
        print(show_unassigned_conversion(master_molID_class_list, step='concatenation')[0][0])

        search_term=="QUIT"
        if search_term == "QUIT":
            search_ok = "OK"
            value = "QUIT"
        search_term = search_term.split(":")
        if (len(search_term) > 4):
            search_ok = "OK"
    return search_term, value


def search_chains(master_molID_class_list, search_term):
    """
    Search the master MolID (entity name) class list for entries that match the provided
    search term, identifying and copying relevant chain IDs and their associated information.

    Parameters
    ----------
    master_molID_class_list : list
        A list of `MolID` class instances representing the master data structure of entity
        names and associated chain IDs.
    search_term : list of str
        A list containing the search criteria. The elements should be in the following order:
        - search_term_file : str
            The name or part of the name of the file to search. An empty string matches any file.
        - search_term_molID : str
            The entity name to search for. An empty string matches any entity name.
        - search_term_oldchID : str
            The old chain ID to search for. An empty string matches any old chain ID.
        - search_term_newchID : str
            The new chain ID to search for. An empty string matches any new chain ID.
        - search_term_concatorder : str
            The concatenation order to search for. An empty string matches any concatenation order.

    Returns
    -------
    found_molID_class_chID_map : dict
        A dictionary mapping copied `MolID` class instances to lists of chain IDs that
        match the search criteria.
    molID_class_been_copied : dict
        A dictionary mapping original `MolID` class instances to their corresponding
        copied versions.
    """
    search_term_file = search_term[0]
    search_term_molID = search_term[1]
    search_term_oldchID = str(search_term[2])
    search_term_newchID = str(search_term[3])
    search_term_concatorder = search_term[4]

    found_molID_class_chID_map = {}
    molID_class_been_copied = {}
    # Perform search
    # Hits to the search have to be copied using deepcopy which will
    # create a distinct object that is an exact copy of another object
    # Modifications will be made to the copy to make sure they do not
    # introduce new errors
    for molID_class in master_molID_class_list:
        if (search_term_file in molID_class.file_name) or (search_term_file == ""):
            for molID in molID_class.molID_chID:
                if (search_term_molID in molID) or (search_term_molID == ""):
                    for chID in molID_class.molID_chID[molID]:
                        if (search_term_oldchID == chID) or (search_term_oldchID == ""):
                            newchID = molID_class.chID_newchID_map[chID]
                            if (search_term_newchID == newchID) or (search_term_newchID == ""):
                                if chID in molID_class.concat_order:
                                    this_concatord = molID_class.concat_order[chID]
                                else:
                                    this_concatord = 0
                                if (search_term_concatorder == str(this_concatord)) or (search_term_concatorder == ""):
                                    print(molID_class.file_name + ":" + molID + ":"
                                          + chID + ":" +
                                          molID_class.chID_newchID_map[chID] + ":"
                                          + str(this_concatord))
                                    if molID_class in molID_class_been_copied:
                                        found_molID_class_chID_map[molID_class_been_copied[molID_class]].append(chID)
                                    else:
                                        test_molID_class = copy.deepcopy(molID_class)
                                        molID_class_been_copied[molID_class] = test_molID_class
                                        found_molID_class_chID_map[test_molID_class] = [chID]
    return found_molID_class_chID_map, molID_class_been_copied

def edit_chain_order(found_map, newchID, action='try'):
    """
    Edits the chain ID or concatenation order for entity name classes based on the specified action.

    This function updates the chain IDs or concatenation order for the entity name classes
    (`molID_class`) in the provided `found_map`. Depending on the action specified, it can
    try a new chain ID, update the concatenation order, or accept the current order as complete.

    Parameters:
    -----------
    found_map : dict
        A dictionary where the keys are `molID_class` objects and the values are lists of
        chain IDs that need to be updated.

    newchID : str
        The new chain ID or concatenation order to be assigned to the chain IDs in `found_map`.

    action : str, optional
        The action to be performed on the chain IDs:
        - 'try' : Attempts to assign `newchID` as the new chain ID and checks for any
                  existing concatenations.
        - 'update' : Updates the concatenation order with the provided `newchID`.
        - 'accept' : Forces the chain ID to be marked as complete.
        Default is 'try'.

    Returns:
    --------
    found_map : dict
        The updated `found_map` with the new chain IDs or concatenation orders applied.
    """
    for molID_class in found_map:
        for chID in found_map[molID_class]:
            if(action=='try'):
                molID_class.chID_newchID_map[chID] = newchID
                #HAVE TO FORCE EXISTING CONCATENATIONS BACK TO OK
                molID_class.check_for_concatenations()
            elif(action=='update'):
                molID_class.update_concat_order(chID, newchID)
            elif(action=='accept'):
                molID_class.force_complete_order(chID, True)
    return found_map

def print_conflicts(found_map):
    """
    Prints any conflicting chain ID assignments that would result from updating to a new chain ID.

    This function checks for any incomplete or conflicting chain ID assignments within the
    provided map of entity name classes (`molID_class`). It identifies and prints details
    of these conflicts, ensuring that each chain ID is only reported once.

    Parameters:
    -----------
    found_map : dict
        A dictionary where the keys are `molID_class` objects and the values are lists of
        chain IDs that are being considered for updates. The function checks these classes
        for conflicts related to incomplete chain ID assignments.

    Returns:
    --------
    None
    """
    print("Updating this new chain ID will lead to the following conflicting assignments")
    was_printed = {}
    for molID_class in found_map:
        for molID in molID_class.molID_chID:
            for chID in molID_class.molID_chID[molID]:
                if molID_class.complete_order[chID] is False:
                    if chID not in was_printed:
                        was_printed[chID] = 1
                        print(molID_class.file_name + ":" + molID + ":" + chID + ":" + molID_class.chID_newchID_map[chID] + ":" + str(molID_class.concat_order[chID]))

def accept_newchain(masterlist, found_map):
    """
    Accepts and updates the master list with new chain assignments.

    This function updates the master list of entity name classes by replacing
    the existing entries with updated versions from  `found_map`.
    It ensures that each updated MolID class is only added once to the
    master list.

    Parameters:
    -----------
    masterlist : list
        A list of `molID_class` objects representing the current state of
        entity names and their chain IDs.

    found_map : list
        A list of updated `molID_class` objects that contain new chain ID
        assignments to be integrated into the `masterlist`.

    Returns:
    --------
    masterlist : list
        The updated list of `molID_class` objects, reflecting the new chain ID
        assignments.
    """
    usage = {}
    for updated_molID_class in found_map:
        for molID_class in masterlist:
            if (molID_class.file_name == updated_molID_class.file_name):
                if updated_molID_class not in usage:
                    usage[updated_molID_class] = 1
                    masterlist.remove(molID_class)
                    masterlist.append(updated_molID_class)
    return masterlist

#################
# FINALIZE STEP #
#################

def masterlist_to_pdb(filelist, masterlist, target_dir=None):
    """
    Updates CIF files with updated chain IDs based on a master list.

    Function processes a list of CIF files, updating the chain IDs according to the mappings provided in the
    `masterlist`. It renumbers residue IDs if required due to chain concatenations and writes the modified lines to new
    CIF(s) in the specified target directory.

    Parameters:
    -----------
    filelist : list
        A list of file paths to CIF files that need to be processed.

    masterlist : list
        A list of molID class objects, each containing chain ID mappings (`chID_newchID_map`) and residue renumbering
        information (`concat_order`).

    target_dir : str, optional
        The directory where the updated CIF files will be saved. If not provided, the files will be saved in the current
        directory.

    Returns:
    --------
    None
    """
    for my_files in filelist:
        newciffilename=target_dir+'/'+my_files.split('/')[-1]
        with open(my_files) as myfile:
            with open(newciffilename, 'w') as newciffile:
                for molID_class in masterlist:
                    if (molID_class.file_name == myfile.name):
                        for line in myfile:
                            if (line[0:4] == "ATOM") or (line[0:6]=="HETATM"):
                                # Chains outside map should not exist but just in case
                                line_split = line.split()
                                if line_split[6] in molID_class.chID_newchID_map: # FAPA: CHANGING LINE_SPLIT[17] TO 6
                                    # Residues have to be renumbered due to concatenations
                                    if line_split[6] in molID_class.concat_order:
                                        residue_offset = (molID_class.concat_order[line_split[6]] - 1) * 50000 # FAPA changed to 50000 from 1000
                                        new_resinum = int(line_split[15]) + int(residue_offset)
                                        heresthenew_resinum = int(new_resinum)
                                        newline = line_split[0] + " " + line_split[1] + " " + line_split[2] + " " + line_split[3] + " " + line_split[4] + " " + line_split[5] + " " + molID_class.chID_newchID_map[line_split[6]] + " " + line_split[7] + " " + line_split[8] + " " + line_split[9] + " " + line_split[10] + " " + line_split[11] + " " + line_split[12] + " " + line_split[13] + " " + line_split[14] + " " + str(heresthenew_resinum) + " " + line_split[16] + " " + molID_class.chID_newchID_map[line_split[6]] + " " + line_split[18] + " " + line_split[19] + "\n"
                                        newciffile.write(newline)
                                    else:
                                        newline = line_split[0] + " " + line_split[1] + " " + line_split[2] + " " + line_split[3] + " " + line_split[4] + " " + line_split[5] + " " + molID_class.chID_newchID_map[line_split[6]] + " " + line_split[7] + " " + line_split[8] + " " + line_split[9] + " " + line_split[10] + " " + line_split[11] + " " + line_split[12] + " " + line_split[13] + " " + line_split[14] + " " + line_split[15] + " " + line_split[16] + " " + molID_class.chID_newchID_map[line_split[6]] + " " + line_split[18] + " " + line_split[19] + "\n"
                                        newciffile.write(newline)

                                else:
                                    newciffile.write(line)

                            else:
                                newciffile.write(line)
