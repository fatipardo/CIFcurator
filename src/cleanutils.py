import os, glob
import re
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

# OCT2026: PDBx/mmCIF dictionary that the output files follow, declared in their _audit_conform category
MMCIF_DICT_NAME = 'mmcif_pdbx.dic'
MMCIF_DICT_VERSION = '5.418'
MMCIF_DICT_LOCATION = 'https://mmcif.wwpdb.org/dictionaries/ascii/mmcif_pdbx_v50.dic'

def process(projdir=None, step='clean', source='raw_bank', target='clean_bank', pdbformat='.cif', verbose=True):
    """
    Processes all CIF files in the source directory through the cleaning step, and then saves the results
    to the target directory. OCT2026: 'clean' and 'simplify' are now a single step (see clean_cif);
    step='simplify' is still accepted and does the same.

    Parameters:
    -----------
    projdir : str, optional
        The project directory containing source and target subdirectories. If none, the current directory is used.
    step : str, optional
        The processing step to apply to each CIF(s). The default, 'clean', rewrites the CIF(s) with only a limited
        set of categories, one file per biological assembly.
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

            # OCT2026: clean and simplify merged into one step
            if(step=='clean' or step=='simplify'):
                clean_cif(input_cif, output_cif, pdbformat)

def clean_cif(oldfile, newfile, pdbformat='.cif'):
    """
    Rewrites a CIF with only a limited set of categories, one file per biological assembly.

    OCT2026: the former 'clean' and 'simplify' steps, merged. The output keeps the categories and key items
    that the data in it require (atom_type, chem_comp, struct_asym, entity_poly, entity_poly_seq,
    pdbx_poly_seq_scheme; citation.id, citation_author.citation_id/ordinal, exptl.entry_id), so that it is
    consistent with the mmCIF dictionary. As before, entity names are written in upper case, without
    colons or single quotes, and _atom_site keeps the 20 columns that Steps 3-5 read.

    A wwPDB assembly file, as downloaded in Step 1 (PDBID+NN.cif, no _pdbx_struct_assembly_gen), already
    holds a single assembly: it is rewritten whole, under the same name. A full PDB entry is split into
    one file per assembly (PDBID+0N.cif, or PDBID+00.cif if it has only one).

    Parameters:
    -----------
    oldfile : str
        Path to the original CIF, as downloaded in Step 1.
    newfile: str
        Path where the new CIF will be saved. Multiple files are created if a full entry has more than one
        biological assembly.
    pdbformat : str, optional
        The file format extension used when saving the new CIF. Default is '.cif'.

    Returns:
    -----------
    None
    """
    mmcif_dict = MMCIF2Dict(oldfile)

    # OCT2026: newfile_stem moved up: it is also the name of a rewritten assembly file
    # OCT2026: strip only a trailing pdbformat. re.sub(pdbformat, ...) treated ".cif" as a regex
    # and removed every "?cif" anywhere in the path (e.g. /tmp/specific_bank -> /tmp/spic_bank)
    newfile_stem = newfile[:-len(pdbformat)] if newfile.endswith(pdbformat) else newfile #OCT2026

    # OCT2026: oldfile is a wwPDB assembly file: written whole, with all atoms. Its _entry.id is XXXX:
    # the PDB ID is taken from the filename instead. Extended IDs (pdb_00001aw1) are lower case,
    # as in the wwPDB beta archive entries; legacy IDs (1AW1) are upper case
    if '_pdbx_struct_assembly_gen.assembly_id' not in mmcif_dict: #OCT2026
        oldfile_stem = oldfile[:-len(pdbformat)] if oldfile.endswith(pdbformat) else oldfile
        pdbid = os.path.basename(oldfile_stem).split('+')[0]
        entry_id = "'"+(pdbid.lower() if pdbid.lower().startswith('pdb_') else pdbid.upper())+"'"
        outputs = [(newfile_stem, None)]
    else:
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


        # Changes the list format to str
        # Update 02/20/2026 : Paulina Removed brackets which would surround ' _entry.id'
        entry_id = str(mmcif_dict['_entry.id']).replace("]","").replace("[","")

        # This section is necessary to print the biological assemblies on separate files
        BioAssembly = mmcif_dict['_pdbx_struct_assembly_gen.asym_id_list']

        # One output per assembly: (file name without extension, label_asym_ids to keep)
        outputs = []
        for assembly in assembly_id_list:
            if (len(assembly_id_list)==1):
                newciffilename = newfile_stem+"+00" #OCT2026
            else:
                newciffilename = newfile_stem+"+0"+str(assembly) #OCT2026
            outputs.append((newciffilename, BioAssembly[int(assembly)-1].split(',')))

    # OCT2026: quotes a value as mmCIF requires (bare, '...', "..." or a ;-delimited text field)
    def cif_value(value):
        if value == '':
            return "''"
        if '\n' in value or ("'" in value and '"' in value and re.search(r"\s", value)):
            return "\n;" + value + "\n;\n"
        if not re.search(r"\s", value) and value[0] not in "_#$'\"[];" and value.lower() not in ('loop_', 'stop_', 'global_') and not value.lower().startswith(('data_', 'save_')):
            return value
        if "'" not in value:
            return "'" + value + "'"
        return '"' + value + '"'

    for newciffilename, asym_ids in outputs:
        with open(newciffilename+pdbformat, 'w') as newciffile:
            # OCT2026: the block name is the file name, not its full path
            newciffile.write("data_"+os.path.basename(newciffilename)+"\n")

            # Writes entry.id
            newciffile.write("#\n")

            entryid = '_entry.id   ' + entry_id
            newciffile.write(entryid + "\n")

            # OCT2026: Writes the dictionary the file follows
            newciffile.write("#\n")
            newciffile.write("_audit_conform.dict_name       " + MMCIF_DICT_NAME + "\n")
            newciffile.write("_audit_conform.dict_version    " + MMCIF_DICT_VERSION + "\n")
            newciffile.write("_audit_conform.dict_location   " + MMCIF_DICT_LOCATION + "\n")

            # Write Audit category
            # OCT2026: with its key items citation_id and ordinal
            newciffile.write("#\n")
            newciffile.write("loop_\n")
            newciffile.write("_citation_author.citation_id\n")
            newciffile.write("_citation_author.name\n")
            newciffile.write("_citation_author.ordinal\n")

            if '_citation_author.name' in mmcif_dict:
                L = mmcif_dict['_citation_author.name']
            else:
                L = ["???"]
            n = len(L)
            L0 = mmcif_dict.get('_citation_author.citation_id', ['primary'] * n)
            L2 = mmcif_dict.get('_citation_author.ordinal', [str(i + 1) for i in range(n)])

            for i in range(n):
                newciffile.write(cif_value(L0[i]) + " '" + re.sub("'", "", L[i]) + "' " + L2[i] + "\n")

            # Writes Citation category
            # OCT2026: with its key item id
            newciffile.write("#" + "\n")
            newciffile.write("loop_" + "\n")
            newciffile.write("_citation.id" + "\n")
            newciffile.write("_citation.title" + "\n")
            newciffile.write("_citation.year" + "\n")
            newciffile.write("_citation.pdbx_database_id_DOI" + "\n")
            L1 = mmcif_dict['_citation.title']
            L2 = mmcif_dict['_citation.year']
            L3 = mmcif_dict['_citation.pdbx_database_id_DOI']
            L0 = mmcif_dict.get('_citation.id', ['primary'] + [str(i + 1) for i in range(1, len(L1))])
            for i in range(len(L1)):
                newciffile.write(cif_value(L0[i]) + " '" + re.sub("\n"," ",re.sub("'", "", L1[i])) + "' " + L2[i] + " " + cif_value(L3[i]) + "\n") #FAPA

            # Writes Experimental method category
            # OCT2026: with its key item entry_id. _exptl.resolution is not a dictionary item and is not written
            # OCT2026: computed structure models (e.g. AlphaFold DB) have no experimental method, and _exptl
            # is not written for them ('THEORETICAL MODEL' is only allowed as a supporting method)
            if '_exptl.method' in mmcif_dict:
                L1 = mmcif_dict['_exptl.method']
            elif '_refine_hist.pdbx_refine_id' in mmcif_dict:
                L1 = mmcif_dict['_refine_hist.pdbx_refine_id']
            elif '_refine.pdbx_refine_id' in mmcif_dict:
                L1 = mmcif_dict['_refine.pdbx_refine_id']
            else:
                L1 = []
            if not isinstance(L1, list):
                L1 = [L1]
            if L1:
                newciffile.write("#" + "\n")
                newciffile.write("loop_" + "\n")
                newciffile.write("_exptl.entry_id" + "\n")
                newciffile.write("_exptl.method" + "\n")
            for i in range(len(L1)):
                newciffile.write(entry_id + " '" + L1[i] + "'\n")

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

            # OCT2026: categories that the atoms refer to. Rows of struct_asym and pdbx_poly_seq_scheme are kept
            # only for the chains of this assembly. Steps 3-5 rename chains and residues in these two
            # categories too, and read their rows by position: keep the column order.
            kept_categories = [('_entity_poly', ['entity_id', 'type', 'nstd_linkage', 'nstd_monomer'], None),
                               ('_entity_poly_seq', ['entity_id', 'num', 'mon_id', 'hetero'], None),
                               ('_chem_comp', ['id', 'type', 'mon_nstd_flag', 'name', 'pdbx_synonyms', 'formula',
                                               'formula_weight'], None),
                               ('_struct_asym', ['id', 'entity_id'], 'id'),
                               ('_pdbx_poly_seq_scheme', ['asym_id', 'entity_id', 'seq_id', 'mon_id', 'ndb_seq_num',
                                                          'pdb_seq_num', 'auth_seq_num', 'pdb_mon_id', 'auth_mon_id',
                                                          'pdb_strand_id', 'pdb_ins_code', 'hetero'], 'asym_id')]
            for category, items, asym_item in kept_categories:
                if category + '.' + items[0] not in mmcif_dict:
                    continue
                n = len(mmcif_dict[category + '.' + items[0]])
                columns = [mmcif_dict.get(category + '.' + item, ['?'] * n) for item in items]
                newciffile.write("#" + "\n")
                newciffile.write("loop_" + "\n")
                for item in items:
                    newciffile.write(category + "." + item + "\n")
                for i in range(n):
                    if asym_ids is None or asym_item is None or mmcif_dict[category + '.' + asym_item][i] in asym_ids:
                        newciffile.write(" ".join(cif_value(column[i]) for column in columns) + "\n")

            # OCT2026: one atom_type row per element in the atoms written below
            newciffile.write("#" + "\n")
            newciffile.write("loop_" + "\n")
            newciffile.write("_atom_type.symbol" + "\n")
            L3 = mmcif_dict['_atom_site.type_symbol']
            L7 = mmcif_dict['_atom_site.label_asym_id']
            for symbol in sorted(set(L3[i] for i in range(len(L3)) if asym_ids is None or L7[i] in asym_ids)):
                newciffile.write(symbol + "\n")

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

            for i in range(len(L1)):
                if asym_ids is None or L7[i] in asym_ids:
                    newciffile.write(L1[i] + " " + L2[i] + " " + L3[i] + ' "' + L4[i] + '" ' + L5[i] + " " + L6[i] + " " + L7[i] + " " + L8[i] + " " + L9[i] + " " + L10[i] + " " + L11[i] + " " + L12[i] + " " + L13[i] + " " + L14[i] + " " + L15[i] + " " + L16[i] + " " + L17[i] + " " + L18[i] + ' "' + L19[i] + '" ' + L20[i] + "\n")
            newciffile.write("#" + "\n")
