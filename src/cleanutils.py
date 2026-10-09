import os, glob
import re
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

def process(projdir=None, step='clean', source='raw_bank', target='clean_bank', pdbformat='.cif', verbose=True):
    """
    Processes all CIF files in the source directory through one of the processing steps,
    and then saves the results to the target directory. The specified steps include, 'clean' and 'simplify'

    Parameters:
    -----------
    projdir : str, optional
        The project directory containing source and target subdirectories. If none, the current directory is used.
    step : str, optional
        The processing step to apply to each CIF(s). If none, the current step is 'clean', which rewrites the CIF(s)
        including only a limited set of data blocks.
    source : str, optional
        The subdirectory within the project directory where the raw CIF(s) are in. The default is the
        subdirectory titled 'raw_bank'
    target : str, optional
        The subdirectory within the project directory where processed CIF(s) will be saved. The default is the
        subdirectory titled 'clean_bank'
    pdbformat : str, optional
        The file extension format for CIF(s). Thr default is '.cif'.
    verbose : bool, optional
        If True, progress is printed to the console. Default is true.

    Returns:
    -----------
    None
    """
    if projdir is not None:
        source_dir = projdir+'/'+source
        target_dir = projdir+'/'+target
        input_list = glob.glob(source_dir+'/*'+pdbformat)
        i=0

        for input_cif in input_list:
            cif_name=os.path.basename(input_cif)

            if verbose:
                i+=1
                print('[{0}/{1}]: {2}'.format(i,len(input_list),cif_name))

            output_cif=target_dir+'/'+cif_name

            if(step=='clean'):
                if os.path.isfile(output_cif):
                    os.remove(output_cif)
                clean_cif(input_cif, output_cif)

            elif(step=='simplify'):
                # missing line: remove all assembly cif already created
                simplify_cif(input_cif, output_cif, pdbformat)

def simplify_cif(oldfile, newfile, pdbformat):
    """
    Separate all biological assemblies in a CIF into separate files.

    Parameters:
    -----------
    oldfile : str
        Path to the original CIF(s) needed to be simplified.
    newfile: str
        Path where the new, simplified CIF(s) will be saved. The function creates multiple files if there are
        more than one biological assemblies.
    pdbformat : str
        The file format extension used when saving the new CIF.

    Returns:
    -----------
    None
    """
    mmcif_dict = MMCIF2Dict(oldfile)

    # Create map from asym_id to assembly_id
    # Convert assembly_id to a list, as it can be either a string or a list

    asym_assembly_map = {}
    assembly_id = mmcif_dict['_pdbx_struct_assembly_gen.assembly_id']

    if not isinstance(assembly_id, list):
        assembly_id_list = []
        asym_id_list = []
        assembly_id_list.append(assembly_id)
        asym_id_list.append(mmcif_dict['_pdbx_struct_assembly_gen.asym_id_list'])

    else:
        assembly_id_list = assembly_id
        asym_id_list = mmcif_dict['_pdbx_struct_assembly_gen.asym_id_list']

    # Convert asym_id entry into a list of asym_ids
    for i in range(len(assembly_id_list)):
        asym_id = asym_id_list[i]
        asym_id = asym_id.strip()
        asym_id = re.sub(' ', '', asym_id)
        asym_id = asym_id.split(',')
        for ident in asym_id:
            asym_assembly_map[ident] = assembly_id_list[i]

    # OCT2026: strip only a trailing pdbformat. re.sub(pdbformat, ...) treated ".cif" as a regex
    # and removed every "?cif" anywhere in the path (e.g. /tmp/specific_bank -> /tmp/spic_bank)
    newfile_stem = newfile[:-len(pdbformat)] if newfile.endswith(pdbformat) else newfile #OCT2026

    for assembly in assembly_id_list:

        if (len(assembly_id_list)==1):
            newciffilename = newfile_stem+"+00" #OCT2026
        else:
            newciffilename = newfile_stem+"+0"+str(assembly) #OCT2026

        newciffile = open(newciffilename+pdbformat, 'w')
        newciffile.write("data_"+newciffilename+"\n")

        # Writes entry.id
        newciffile.write("#\n")

        # Changes the list format to str
        # Update 02/20/2026 : Paulina Removed brackets which would surround ' _entry.id'
        L = str(mmcif_dict['_entry.id']).replace("]","").replace("[","")
        entryid = '_entry.id   ' + L
        newciffile.write(entryid + "\n")

        # Write Audit category
        newciffile.write("#\n")
        newciffile.write("loop_\n")
        newciffile.write("_citation_author.name\n")

        if '_citation_author.name' in mmcif_dict:
            L = mmcif_dict['_citation_author.name']
        else:
            L = "???"


        if isinstance(L, list):
            for i in L:
                newciffile.write("'" + re.sub("'", "", i) + "'" + "\n")
        else:
            newciffile.write("'" + re.sub("'", "", L) + "'" + "\n")

        # Writes Citation category
        newciffile.write("#" + "\n")
        newciffile.write("loop_" + "\n")
        newciffile.write("_citation.title" + "\n")
        newciffile.write("_citation.year" + "\n")
        newciffile.write("_citation.pdbx_database_id_DOI" + "\n")
        L1 = mmcif_dict['_citation.title']
        L2 = mmcif_dict['_citation.year']
        L3 = mmcif_dict['_citation.pdbx_database_id_DOI']
        if isinstance(L1, list):
            for i in range(len(L1)):
                newciffile.write("'" + re.sub("\n"," ",re.sub("'", "", L1[i])) + "' " + L2[i] + " " + L3[i] + "\n") #FAPA
        else:
            newciffile.write("'" + re.sub("\n"," ",re.sub("'", "", L1[i])) + "' " + L2[i] + " " + L3[i] + "\n") #FAPA

        # Writes Resolution category
        newciffile.write("#" + "\n")
        newciffile.write("loop_" + "\n")
        newciffile.write("_exptl.method" + "\n")
        newciffile.write("_exptl.resolution" + "\n")
        if '_exptl.method' in mmcif_dict:
            L1 = mmcif_dict['_exptl.method']
        elif '_refine_hist.pdbx_refine_id' in mmcif_dict:
            L1 = mmcif_dict['_refine_hist.pdbx_refine_id']
        else:
            L1 = mmcif_dict['_refine.pdbx_refine_id']
        if '_refine.ls_d_res_high' in mmcif_dict:
            L2 = mmcif_dict['_refine.ls_d_res_high']
        elif '_em_3d_reconstruction.resolution' in mmcif_dict:
            L2 = mmcif_dict['_em_3d_reconstruction.resolution']
        elif '_refine_hist.d_res_high' in mmcif_dict:
            L2 = mmcif_dict['_refine_hist.d_res_high']
        else:
            L2 = '????'
        if isinstance(L1, list) and isinstance(L2, list):
            for i in range(len(L1)):
                newciffile.write("'" + L1[i] + "' " + L2[i] + " " + "\n")
        elif isinstance(L1, list) and not isinstance(L2,list):
            newciffile.write("'" + L1[0] + "' " + L2 + " " + "\n")
        elif not isinstance(L1,list) and isinstance(L2,list):
            newciffile.write("'" + L1 + "' " + L2[0] + " " + "\n")
        else:
            newciffile.write("'" + L1 + "' " + L2 + " " + "\n")

        # Writes Entity category
        # PAULINA 05/08/26
        newciffile.write("#" + "\n")
        newciffile.write("loop_" + "\n")
        newciffile.write("_entity.id" + "\n")
        newciffile.write("_entity.type" + "\n")
        newciffile.write("_entity.src_method" + "\n")
        newciffile.write("_entity.pdbx_description" + "\n")
        L1 = mmcif_dict['_entity.id']
        n = len(L1)
        L2 = mmcif_dict.get('_entity.type', ['polymer'] * n)
        L3 = mmcif_dict.get('_entity.src_method', ['?'] * n)
        L4 = mmcif_dict['_entity.pdbx_description']
        for i in range(n):
            L4[i] = L4[i].upper()
            L4[i] = L4[i].replace(":", "")
            newciffile.write(L1[i] + " " + L2[i] + " " + L3[i] + " '" + L4[i].replace("'", "") + "'\n")

        # Writes the coordinate portion of the file
        newciffile.write("#" + "\n")
        newciffile.write("loop_" + "\n")
        newciffile.write("_atom_site.group_PDB" + "\n")
        newciffile.write("_atom_site.id" + "\n")
        newciffile.write("_atom_site.type_symbol" + "\n")
        newciffile.write("_atom_site.label_atom_id" + "\n")
        newciffile.write("_atom_site.label_alt_id" + "\n")
        newciffile.write("_atom_site.label_comp_id" + "\n")
        newciffile.write("_atom_site.label_asym_id" + "\n")
        newciffile.write("_atom_site.label_entity_id" + "\n")
        newciffile.write("_atom_site.label_seq_id" + "\n")
        newciffile.write("_atom_site.pdbx_PDB_ins_code" + "\n")
        newciffile.write("_atom_site.Cartn_x" + "\n")
        newciffile.write("_atom_site.Cartn_y" + "\n")
        newciffile.write("_atom_site.Cartn_z" + "\n")
        newciffile.write("_atom_site.occupancy" + "\n")
        newciffile.write("_atom_site.B_iso_or_equiv" + "\n")
        newciffile.write("_atom_site.auth_seq_id" + "\n")
        newciffile.write("_atom_site.auth_comp_id" + "\n")
        newciffile.write("_atom_site.auth_asym_id" + "\n")
        newciffile.write("_atom_site.auth_atom_id" + "\n")
        newciffile.write("_atom_site.pdbx_PDB_model_num" + "\n")
        L1 = mmcif_dict['_atom_site.group_PDB']
        L2 = mmcif_dict['_atom_site.id']
        L3 = mmcif_dict['_atom_site.type_symbol']
        L4 = mmcif_dict['_atom_site.label_atom_id']
        L5 = mmcif_dict['_atom_site.label_alt_id']
        L6 = mmcif_dict['_atom_site.label_comp_id']
        L7 = mmcif_dict['_atom_site.label_asym_id']
        L8 = mmcif_dict['_atom_site.label_entity_id']
        L9 = mmcif_dict['_atom_site.label_seq_id']
        L10 = mmcif_dict['_atom_site.pdbx_PDB_ins_code']
        L11 = mmcif_dict['_atom_site.Cartn_x']
        L12 = mmcif_dict['_atom_site.Cartn_y']
        L13 = mmcif_dict['_atom_site.Cartn_z']
        L14 = mmcif_dict['_atom_site.occupancy']
        L15 = mmcif_dict['_atom_site.B_iso_or_equiv']
        L16 = mmcif_dict['_atom_site.auth_seq_id']
        L17 = mmcif_dict['_atom_site.auth_comp_id']
        L18 = mmcif_dict['_atom_site.auth_asym_id']
        L19 = mmcif_dict['_atom_site.auth_atom_id']
        L20 = mmcif_dict['_atom_site.pdbx_PDB_model_num']

        # This section is necessary to print the biological assemblies on separate files
        BioAssembly = mmcif_dict['_pdbx_struct_assembly_gen.asym_id_list']

        for i in range(len(L1)):
            if (L7[i] in BioAssembly[int(assembly)-1].split(',')):
                newciffile.write(L1[i] + " " + L2[i] + " " + L3[i] + ' "' + L4[i] + '" ' + L5[i] + " " + L6[i] + " " + L7[i] + " " + L8[i] + " " + L9[i] + " " + L10[i] + " " + L11[i] + " " + L12[i] + " " + L13[i] + " " + L14[i] + " " + L15[i] + " " + L16[i] + " " + L17[i] + " " + L18[i] + ' "' + L19[i] + '" ' + L20[i] + "\n")
        newciffile.write("#" + "\n")


#
def clean_cif(oldfile, newfile):
    """
    Rewrites CIF, including only a limited set of data blocks.

    Parameters:
    -----------
    oldfile : str
        The path to the original CIF(s) needed to be cleaned.
    newfile : str
        The path where the cleaned CIF(s) will be written.

    Returns:
    -----------
    None
    """
    entry_list = ['_entry.id',
                  '_atom_site.group_PDB',
                  '_citation_author.name',
                  '_citation.title',
                  '_pdbx_struct_assembly_gen.assembly_id',
                  '_entity.pdbx_description',
                  '_exptl.method',
                  '_em_3d_reconstruction.resolution',
                  '_refine_hist.pdbx_refine_id',
                  '_refine.pdbx_refine_id']
    keylength_list = [ 9,
                       20,
                       21,
                       15,
                       37,
                       24,
                       13,
                       32,
                       27,
                       22]
    with open(oldfile) as old_file:
        alllines = []
        linecount = 0
        poundline = 0
        flag = 0
        for line in old_file:
            alllines.append(line)
            if linecount == 0:
                with open(newfile, 'a') as new_file:
                    new_file.write(alllines[0])
            for entry, keylength in zip(entry_list, keylength_list):
                flag = check_and_write_entry(entry, line, alllines, line[0:keylength], flag, range(poundline, linecount), newfile)
            if '#' in line[0]:
                poundline = linecount
            linecount += 1
        with open(newfile, 'a') as new_file:
            new_file.write('#\n')
#
def check_and_write_entry(entry, line, alllines, key, flag, linerange, newfile):
    """
    Checks if a specific entry is present in the current line of a CIF and writes relevant lines to a new file.

    Parameters:
    -----------
    entry : str
        The specific CIF entry to look for in the line (e.g., '_entry.id').
    line : str
        The current line being read from the CIF.
    alllines : list
        A list of all lines read so far from the CIF.
    key : str
        The substring of the current line that is compared to the entry.
    flag : int
        A flag indicating whether the desired entry has been found (1 if found, 0 otherwise).
    linerange : range
        The range of lines from `alllines` to write to the new file if the entry is found.
    newfile : str
        The path to the new CIF where the relevant lines will be written.

    Returns:
    -----------
    flag : int
        The new flag value indicating whether the entry was found or not.
    """
    if entry in key:
        flag = 1
    elif (flag==1) and '#' in line[0]:
        with open(newfile, 'a') as new_file:
            for i in linerange:
                new_file.write(alllines[i])
        flag=0
    return flag
