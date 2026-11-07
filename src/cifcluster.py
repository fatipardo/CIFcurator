import os
import shutil
import subprocess
import pandas as pd
import csv
import string

# ------------------ CLUSTERING ---------------------
def mapping_sequences_count(file_path):
    """
    Analyze sequence and entity frequencies in a FASTA file.

    This function reads sequence and entity annotations from a FASTA file using
    `get_attributes_seqres`, then computes several frequency-based summaries:

    - The most frequent sequence in the file.
    - The most frequent entity associated with that top sequence.
    - All unique sequences, sorted from most to least frequent.
    - All unique entities, sorted from most to least frequent.
    - The frequency count of the top sequence.

    Parameters
    ----------
    file_path(str): Path to the FASTA file

    Returns
    -------
    dict
        A dictionary containing:

        - **"UNIQUE_ENTITIES_SORTED"** : list of str
          All unique entity labels sorted by descending frequency.

        - **"UNIQUE_SEQUENCES_SORTED"** : list of str
          All unique sequences sorted by descending frequency.

        - **"TOP_ENTITY"** : str
          The entity most frequently associated with the top sequence.

        - **"TOP_SEQUENCE"** : str
          The most frequent sequence in the file.

        - **"FREQUENCY_OF_TOP_SEQUENCE"** : int
          Number of occurrences of the top sequence.

        - **"ENTITY_TO_SEQUENCES"** : dict
          Mapping of unique entities to their corresponding sets of sequences.

        - **"PDB_ENTRIES"** : list
          List of PDB chain IDs extracted from the file.

        - **"PDB_ENTRY_TO_ENTITY"** : dict
          Mapping from PDB chain IDs to entity labels.

        - **"PDB_ENTRY_TO_SEQUENCE"** : dict
          Mapping from PDB chain IDs to sequences.

    """
    # Returns of list of sequences and entities from FASTA file
    sequences = get_attributes_seqres(file_path, attribute='SEQUENCE')
    entities = get_attributes_seqres(file_path, attribute='ENTITY_LABEL')
    pdbchid_entity = get_attributes_seqres(file_path, attribute='PDBID_CHID')
    
    sequence_counts = {}
    entity_counts = {}
    entity_to_sequences = {}
    # creates dictionary where = sequence/entity : occurrence_count
    for seq, ent in zip(sequences, entities):
        entity_to_sequences.setdefault(ent, set()).add(seq)
        sequence_counts[seq] = sequence_counts.get(seq, 0) + 1
        entity_counts[ent] = entity_counts.get(ent, 0) + 1

    # I order the dictionaries to have the most frequent to least frequent entity/sequence first
    temp_sorted_unique_entities = sorted(entity_counts.items(), key=lambda x: x[1], reverse=True)
    sorted_unique_entities = [ent[0] for ent in temp_sorted_unique_entities]

    temp_sorted_unique_sequences = sorted(sequence_counts.items(), key=lambda x: x[1], reverse=True)
    sorted_unique_sequences = [seq[0] for seq in temp_sorted_unique_sequences]

    most_freq_seq = max(sequence_counts, key=sequence_counts.get)
    frequency_of_top_seq = sequence_counts[most_freq_seq]

    # I create dictionary mapping the top entity to its corresponding sequence
    top_ent_4_top_seq = {}
    for seq, ent in zip(sequences, entities):
        if seq == most_freq_seq:
            top_ent_4_top_seq[ent] = top_ent_4_top_seq.get(ent, 0) + 1
            
    most_freq_entity = max(top_ent_4_top_seq, key=top_ent_4_top_seq.get)
    info_dict = {"UNIQUE_ENTITIES_SORTED":sorted_unique_entities,
                 "UNIQUE_SEQUENCES_SORTED":sorted_unique_sequences,
                 "TOP_ENTITY":most_freq_entity,
                 "TOP_SEQUENCE":most_freq_seq,
                 "FREQUENCY_OF_TOP_SEQUENCE":frequency_of_top_seq,
                 "ENTITY_TO_SEQUENCES":entity_to_sequences,
                 "PDB_ENTRIES":pdbchid_entity,
                 "PDB_ENTRY_TO_ENTITY": dict(zip(pdbchid_entity, entities)),
                 "PDB_ENTRY_TO_SEQUENCE": dict(zip(pdbchid_entity, sequences))}

    return info_dict

def get_attributes_seqres(file_path, attribute='PDBID_CHID', cif=None):
    """
    Extracts specified attributes (e.g., PDB ID, chain ID, sequence, or entity label)
    from a SEQRES FASTA-like file.

    Parameters
    ----------
    file_path(str): Path to the input seqres.txt file.

    attribute(str, optional): Metadata attribute to extract.
        Options are:
            - 'ENTITY_LABEL' : the entity name
            - 'MOL'          : the molecular identifier
            - 'PDBID'        : the PDB ID only
            - 'CHID'         : the chain ID only
            - 'PDBID_CHID'   : PDB ID and chain ID (default)
            - 'HEADER'       : returns meta-data line
            - 'SEQUENCE'     : the sequence

    cif(str, optional): A specific PDB ID to filter results. If provided, only entries with this PDB ID are included.

    Returns
    -------
    attribute_list(list): A list of extracted attributes or sequences, depending on the `attribute` specified.
    """
    attribute_list = []
    attribute = attribute.upper()
    pdb_id_filter = cif.upper()[:12] if cif else None

    # Had to add the encoding parameter due to an error that kept occuring when reading the file.
    with open(file_path, 'r', encoding='latin-1') as seqres:
        for i, line in enumerate(seqres):

            if i % 2 == 0:
                if attribute == "HEADER":
                    attribute_list.append(line.strip())

                meta_data_parts = line.split()
                pdb_id = meta_data_parts[0].split('_') # >PDB_0000100D_A -> [>PDB, 0000100D, A]
                pdb_id = (pdb_id[0][1:] + "_" + pdb_id[1]).upper()# PDB_0000100D
                mol_id = meta_data_parts[1].strip()
                chid = meta_data_parts[0].split('_')[2]
                pdb_chid = pdb_id + '_' + chid
                entity_name = ' '.join(meta_data_parts[3:]).replace("'", "")

                if pdb_id_filter and pdb_id != pdb_id_filter:
                    continue

                elif attribute == 'ENTITY_LABEL':
                    attribute_list.append(entity_name)

                elif attribute == 'MOL':
                    attribute_list.append(mol_id)

                elif attribute == 'CHID':
                    attribute_list.append(chid)

                elif attribute == 'PDBID_CHID':
                    attribute_list.append(pdb_chid)

                elif attribute == 'PDBID':
                    attribute_list.append(pdb_id)

            elif i % 2 == 1:
                sequence = line.strip().upper()

                if pdb_id_filter and pdb_id != pdb_id_filter:
                    continue

                if attribute == 'SEQUENCE':
                    attribute_list.append(sequence)

    return attribute_list

def entry_dict_seqres(seqres_file):
    """
    Builds a dictionary from SEQRES, mapping each PDB ID and chain ID to its corresponding metadata line and sequence.

    Parameters
    ----------
    seqres_file(str): Path to the input seqres.txt file.

    Returns
    -------
    entry_dict(dict): A dictionary where:
            - Keys are strings in the format "PDBID_CHID" (e.g., "1ABC_A").
            - Values are lists containing: [metadata_line (str), sequence_line (str)].

    Notes
    -----
    - All metadata components except the chain ID are capitalized for consistency.
    """
    entry_dict = {}
    with open(seqres_file, 'r') as file:
        for meta_line, sequence_line in zip(file, file):
            # Here we capitalize everything in pdb_chid except for the chain id
            pdbchid_label_fasta = meta_line.split()[0][1:]

            pdb_id = pdbchid_label_fasta[:12].upper()
            ch_id = pdbchid_label_fasta.split('_')[2]
            pdb_chid = pdb_id + '_' + ch_id

            # We obtain all other components in our metaline (e.g mol, length, entity name)
            rest_line_cap = [word.upper() for word in meta_line.split()[1:]]

            # We construct our metaline with our capitalized components
            meta_line = '>' + pdb_chid + ' ' + ' '.join(rest_line_cap) + '\n'

            # We make our dictionary from our capitalzed seqres components
            entry_dict[pdb_chid] = [meta_line, sequence_line.upper()]

    return entry_dict

def map_unique(list1, list2):
    """
    Creates a mapping from each unique element in list2 to a set of corresponding elements from list1.

    Parameters
    -----------
    list1(list): A list of values to be grouped.
    list2(list): A list of keys which will be used for grouping.

    Returns
    --------
    values_map(dict) or None: Dictionary with list2 elements as keys, and list1 elements as values,
                              or None if list lengths do not match.
    """
    zip_values = zip(list1, list2)
    values_map = {}

    if len(list1) != len(list2):
        return print('Both lists must be of equal lengths.')

    else:
        for value1, value2 in zip_values:
            if value2 not in values_map:
                values_map[value2] = set()
            values_map[value2].add(value1)

    return values_map

def cleaning_clusters(clusteringStorage):
    """
    Goes through the clusters.tsv file and gets unique instances of cluster representatives, and maps them to
    their cluster members.

    Parameters
    -----------
    clusteringStorage(str): path to directory containing clusters.tsv file.

    Returns
    --------
    cluster_dict(dict): dictionary with unique cluster representatives as keys and grouped members as values.
    """
    # fix indexing, extract whole thing
    with open(f"{clusteringStorage}/clusters.tsv", 'r') as file:
        label_list = []
        representatives = []
        members = []
        for i, line in enumerate(file):
            label_list.append(line)
        for label in label_list:
            representatives.append(label.split()[0])
            members.append(label.split()[1])

        cluster_dict = map_unique(members, representatives)
    return cluster_dict

def writing_cluster_FASTA(seqres_dict, cluster_dict, clusterFastaStoragePath):
    """
    Created FASTA files out of the clusters made using MMSeqs2.

    Parameters
    ----------
    seqres_dict(dict): dictionary made out of our seqres file utilizing entry_dict_seqres function

    cluster_dict(dict): dictionary built off of the clusters.tsv having used the cleaning_clusters function

    clusterFastaStoragePath(str): directory where outputted files will be stored

    Returns
    -------
    None
    """
    cluster_num = 1
    for key in cluster_dict:
        with open(f'{clusterFastaStoragePath}/{cluster_num}_{key}_cluster.txt', 'w') as file:
            entries_list = list(cluster_dict[key])
            for entry in entries_list:
                lines_to_write = seqres_dict[entry]
                if lines_to_write:
                    file.writelines(lines_to_write)
                else:
                    print(f'{entry} not in seqres.txt')
        cluster_num += 1

def writing_entity_FASTA(clustering_directory, seqres_dict, PDB_ID_list):
    """
    Writes a FASTA file with entries from seqres_dict related to the PDB IDs of interest.

    NOTE: the txt file containing PDB IDs must be in the main directory we are working from

    Parameters
    ----------
    clustering_directory(str): The main clustering directory where sub_seqres.txt will be stored.
    seqres_dict(dict): A dictionary produced by 'entry_dict_seqres()' where:
            - Keys are strings in the format "PDBID_CHID" (e.g., "1ABC_A").
            - Values are lists containing: [metadata_line (str), sequence_line (str)].
    PDB_ID_list(list): List of PDB IDs to filter and write out.

    Returns
    -------
    None
    """
    # seqres_dict dictionary of the big seqres
    matching_labels = [label for label in seqres_dict if label[:12] in PDB_ID_list] # 0000100D_A
    entityFasta = f'sub_seqres.txt'
    with open(f'{clustering_directory}/{entityFasta}', 'w') as file:
        for label in matching_labels:
            lines_to_write = seqres_dict[label]
            if lines_to_write:
                file.writelines(lines_to_write)
            else: print(f'{label} not in seqres.txt')

def writing_labels_file(clusters_directory_path, output_directory):
    """
    Creates a text file where each line contains the name of a cluster file, the most frequent
    entity name for that cluster (along with their counts), and the unique entity names within
    that cluster (also with its count).

    Parameters
    ----------
    clusters_directory_path(str): path to directory containing the cluster FASTA files

    output_directory(str): directory where the newly created txt file will be stored

    Returns
    -------
    None
    """

    def obtaining_label_counts(file_path):
        """
        Returns a dictionary, with the keys being the unique entity names in a FASTA file
        and the values being the counts of those keys in the file.

        Parameters
        -----------
        file_path : path to a FASTA file

        Returns
        --------
        Dictionary with keys being entity name and values being its number of instances
        """
        label_list = get_attributes_seqres(file_path, 'ENTITY_LABEL')
        unique_label_list = set(label_list)
        counts = []

        for unique_label in unique_label_list:
            count = 0
            for label in label_list:
                if label == unique_label:
                    count += 1
            counts.append(count)
        return dict(zip(unique_label_list, counts))

    def obtain_high_label(label_dict):
        """
        Takes dictionary from obtaining_label_counts and returns a dictionary containing the key and value of
        the most present entity label in the FASTA file

        Parameters
        ----------
        label_dict(dict): Dictionary produced by `obtaining_label_counts()`, where keys are entity
        labels and values are their occurrence counts.

        Returns
        -------
        dict: A dictionary containing a single key-value pair: {most_frequent_label: count}

              If the input dictionary contains one or zero entries, an empty dictionary is returned.
        """
        if len(label_dict.keys()) > 1:
            label_dict = dict(sorted(label_dict.items(), key=lambda item: item[1], reverse=True))
            return {list(label_dict.keys())[0]: label_dict[list(label_dict.keys())[0]]}
        else:
            return {}

    lines_to_write = []
    structure_count = 0

    # This chunk removes a file we dont need.
    ds_store_path = os.path.join(clusters_directory_path, '.DS_Store')
    if os.path.exists(ds_store_path):
        os.remove(ds_store_path)

    for file in os.listdir(clusters_directory_path):

        full_path = os.path.join(clusters_directory_path, file)

        # This line counts how main chains are in the current file its reading,
        # and adds up that total to its count. It repeats that for every file in the cluster
        structure_count += len(get_attributes_seqres(full_path, 'chid'))
        test_dict = obtaining_label_counts(full_path)
        high_label = obtain_high_label(test_dict)

        if high_label:
            line = [f'{file}', f'{high_label}', f'{test_dict}']
        else:
            line = [f'{file}', f'{test_dict}']

        lines_to_write.append(line)

    with open(f'{output_directory}/unique_label_names.txt', 'w') as curr_file:
        curr_file.writelines(f'TOTAL CHIANS : {structure_count}\n\n')
        for line in lines_to_write:
            curr_file.writelines(line)
            curr_file.writelines('\n')


def clustering(clusteringStorage, seqres_file):
    """
    Runs MMseqs2 database creation, clustering, and TSV conversion, then cleans up temporary files.

    Parameters
    ----------
    clusteringStorage(str): Path to the directory where clustering outputs and temporary files are managed.
    seqres_file(str): Path to the input SEQRES sequence file.

    Returns
    -------
    None
    """

    def _cleanup_files(prefixes: list[str], tmp_dir: str) -> None:
        """Helper to remove MMseqs database files and temporary directory."""
        # Remove directory
        if os.path.exists(tmp_dir):
            shutil.rmtree(tmp_dir)

        # Remove generated database files matching prefixes
        current_dir = os.listdir(".")
        for filename in current_dir:
            for prefix in prefixes:
                if filename.startswith(prefix):
                    try:
                        os.remove(filename)
                    except OSError:
                        pass

    seqres_file = os.path.abspath(seqres_file)
    seqres_dict = entry_dict_seqres(seqres_file)

    original_dir = os.getcwd()
    db_prefix = "pdbDB"
    cluster_prefix = "clusterDB"
    tmp_dir = "tmp"
    output_tsv = "clusters.tsv"

    output_directory = f"{clusteringStorage}"

    print("cwd:", os.getcwd())
    print("entity_seqres:", seqres_file)
    print("absolute:", os.path.abspath(seqres_file))
    print("exists:", os.path.exists(seqres_file))

    try:
        os.chdir(output_directory)

        # 1. Create DB
        print("Creating MMseqs2 database...")
        subprocess.run(
            ["mmseqs", "createdb", seqres_file, db_prefix], check=True
        )

        # 2. Cluster
        print("Clustering sequences...")
        subprocess.run(
            [
                "mmseqs",
                "cluster",
                db_prefix,
                cluster_prefix,
                tmp_dir,
                "--min-seq-id",
                "0.30",
                "--cluster-mode",
                "2",
                "--cov-mode",
                "1",
                "-c",
                "0.8",
            ],
            check=True,
        )

        # 3. Create TSV
        print("Generating TSV output...")
        subprocess.run(
            [
                "mmseqs",
                "createtsv",
                db_prefix,
                db_prefix,
                cluster_prefix,
                output_tsv,
            ],
            check=True,
        )

    finally:
        # 4. Clean up temporary files safely
        print("Cleaning up intermediate files...")
        _cleanup_files([db_prefix, cluster_prefix], tmp_dir)
        os.chdir(original_dir)

    print(f"Clustering complete! Output saved to {output_tsv}")

    print("Cleaning up data...")

    clusterFastaStoragePath = f'{clusteringStorage}cluster_files/'

    cluster_dict = cleaning_clusters(clusteringStorage)
    writing_cluster_FASTA(seqres_dict, cluster_dict, clusterFastaStoragePath)
    writing_labels_file(clusterFastaStoragePath, clusteringStorage)
    cluster_results_to_csv(clusteringStorage)

def cluster_results_to_csv(clusteringStorage):
    """
    Organizes the data generated from clusters and 'mapping_sequences_count()' into a readable CSV format

    Parameters
    ----------
    clusteringStorage(str): Path to the main directory hosting the cluster files and output results.

    Returns
    -------
    None
    """
    clusterFastaStoragePath = f'{clusteringStorage}/cluster_files/'
    cluster_files = os.listdir(clusterFastaStoragePath)

    col1 = []  # Cluster
    col2 = []  # Top Macromolecule
    col3 = []  # Top Sequence
    col4 = []  # Frequency of that sequence
    col5 = []  # Unique entity names of that cluster
    col6 = []  # unique sequences of that cluster
    col7 = []  # entities to their sequences

    for file in cluster_files:
        file_path = f"{clusterFastaStoragePath}{file}"
        cluster_attributes = mapping_sequences_count(file_path)

        col1.append(file)
        col2.append(cluster_attributes["TOP_ENTITY"])
        col3.append(cluster_attributes["TOP_SEQUENCE"])
        col4.append(cluster_attributes["FREQUENCY_OF_TOP_SEQUENCE"])
        col5.append(cluster_attributes["UNIQUE_ENTITIES_SORTED"])
        col6.append(cluster_attributes["UNIQUE_SEQUENCES_SORTED"])
        col7.append(cluster_attributes["ENTITY_TO_SEQUENCES"])

    cluster_df = pd.DataFrame({"Cluster": col1,
                                "Top_Entity": col2,
                                "Top_Sequence": col3,
                                "Sequence_Frequency": col4,
                                "Unique_Entities": col5,
                                "Unique_Sequences": col6,
                               "Entity_to_sequence": col7})

    cluster_df.to_csv(f"{clusteringStorage}/cluster_info(80).csv", index=False)

# WE should be clustering the RNAs and DNAs seperately from the PROTEINS
def processing_RNA_PROTEIN_cluster(clustering_directory, seqres_file):
    """
    Separates protein and RNA/DNA sequences from a SEQRES file and performs clustering for each group independently.

    Parameters
    ----------
    clustering_directory (str): Path to the directory where clustering output directories are stored.
    seqres_file (str): Path to the input SEQRES file containing PDB IDs, molecule types, metadata headers, and sequences.

    Returns
    -------
    None

    Notes
    -----
    - Sequences are classified as RNA/DNA if their molecule type is "MOL:NA"; all other sequences are classified as proteins.
    - Protein and RNA/DNA sequences are processed and clustered separately.
    - Creates a separate clustering storage directory and cluster FASTA directory for each non-empty group.
    - Generates a PDB list file and a SEQRES file for each group if the corresponding files do not already exist.
    - Skips clustering for groups that contain no sequences.
    - Existing PDB list and SEQRES files are preserved and are not overwritten.
    """
    pdbIDs = get_attributes_seqres(seqres_file, attribute="PDBID")
    MOL = get_attributes_seqres(seqres_file, attribute="MOL")
    SEQUENCE = get_attributes_seqres(seqres_file, attribute="SEQUENCE")
    meta_data_line = get_attributes_seqres(seqres_file, attribute="HEADER")

    rna_pdbs = []
    rna_head = []
    rna_seq = []

    protein_pdbs = []
    protein_head = []
    protein_seq = []

    # Separate RNA/DNA and protein sequences
    for pdbid, header, mol, sequence in zip(
        pdbIDs, meta_data_line, MOL, SEQUENCE
    ):
        if mol == "MOL:NA":
            rna_pdbs.append(pdbid)
            rna_head.append(header)
            rna_seq.append(sequence)
        else:
            protein_pdbs.append(pdbid)
            protein_head.append(header)
            protein_seq.append(sequence)

    groups = [
        ("PROTEIN", protein_pdbs, protein_head, protein_seq),
        ("RNA_DNA", rna_pdbs, rna_head, rna_seq),
    ]

    for org, pdbs, heads, seqs in groups:
        # Skip empty groups
        if not seqs:
            print(f"Skipping {org}: no sequences found.")
            continue

        print(f"\nProcessing {org}: {len(seqs)} sequences")

        clusteringStorage = f"{clustering_directory}/{org}_clustering_storage/"
        os.makedirs(clusteringStorage, exist_ok=True)

        clusterFastaStoragePath = f"{clusteringStorage}cluster_files/"
        os.makedirs(clusterFastaStoragePath, exist_ok=True)

        # Write sequence file if it doesn't exist
        seqres_path = os.path.join(
            clusteringStorage, f"{org}_seqres.txt"
        )

        if not os.path.exists(seqres_path):
            with open(seqres_path, "w") as f:
                for head, seq in zip(heads, seqs):
                    f.write(f"{head}\n{seq}\n")

        # Run clustering only for non-empty groups
        print(f"Cluster Stage: {org}")

        clustering(clusteringStorage, seqres_path)

    print("\nRNA/DNA and protein clustering stages completed.")


def finding_problematic_cluster_results(clustering_directory):
    """
    Analyzes clustering results to identify duplicate entities or chains appearing across multiple clusters,
    generating warning CSV tables and returning mapping dictionaries.

    Parameters
    ----------
    clustering_directory (str): Path to the directory containing clustering storage results.

    Returns
    -------
    tuple of (dict, dict)
        - cluster_pdbchids: Mapping of cluster representatives to their members.
        - outlier_case_entries: Mapping of PDB entries or outlier entities to associated issues/entities.
    """
    pdbchid = get_attributes_seqres(f"{clustering_directory}/sub_seqres.txt", attribute='PDBID_CHID')
    entity_names = get_attributes_seqres(f"{clustering_directory}/sub_seqres.txt", attribute='ENTITY_LABEL')

    pdb_entry_to_entity = dict(zip(pdbchid, entity_names))

    directory_names = [f.name for f in os.scandir(clustering_directory) if f.is_dir()]
    cluster_pdbchids = {}
    pdb_entry_to_files = {}
    entity_to_files = {}

    # ==========================================
    # 3. DIRECTORY & CLUSTER FILE PROCESSING
    # ==========================================
    for directory in directory_names:
        cluster_rep = []
        cluster_mem = []
        cluster_file_path = f"{clustering_directory}/{directory}/clusters.tsv"

        if os.path.exists(cluster_file_path):
            with open(cluster_file_path, "r") as f:
                reader = csv.reader(f, delimiter="\t")
                for row in reader:
                    if len(row) >= 2:
                        cluster_rep.append(row[0])
                        cluster_mem.append(row[1])

        cluster_pdbchids.update(
            map_unique(cluster_mem, cluster_rep))  # Its supposed to be cluster_rep : {mem1, mem2, mem3}

        cluster_files_storage = f"{clustering_directory}/{directory}/cluster_files/"
        if os.path.exists(cluster_files_storage):
            for cluster_file in os.listdir(cluster_files_storage):
                full_file_path = f"{cluster_files_storage}/{cluster_file}"
                info = mapping_sequences_count(full_file_path)

                for pdb_entry in info.get("PDB_ENTRIES", []):
                    pdb_entry_to_files.setdefault(pdb_entry, []).append(full_file_path)

                for entity in info.get("UNIQUE_ENTITIES_SORTED", []):
                    entity_to_files.setdefault(entity, []).append(full_file_path)

    pdb_entries_in_multiple_files = {
        k: v for k, v in pdb_entry_to_files.items() if len(v) > 1
    }

    entities_in_multiple_files = {
        k: v for k, v in entity_to_files.items() if len(v) > 1
    }

    # Maps pdb_id_chid -> set of associated entity names
    # Using sets prevents entity names from being overwritten
    outlier_case_entries = {}

    warning_message_pdb_entries_multiple_files = {}

    # Identify entities appearing in multiple cluster files
    if entities_in_multiple_files:
        for key, lst in entities_in_multiple_files.items():
            for cluster_file in lst:

                for entity, pdb_id in zip(
                        get_attributes_seqres(cluster_file, attribute='ENTITY_LABEL'),
                        get_attributes_seqres(cluster_file, attribute='PDBID_CHID')
                ):
                    if entity == key:
                        outlier_case_entries.setdefault(pdb_id, set()).add(entity)

    # Identify PDB chain IDs appearing in multiple cluster files
    if pdb_entries_in_multiple_files:
        for key in pdb_entries_in_multiple_files:
            outlier_case_entries.setdefault(key, set()).add(
                pdb_entry_to_entity[key]
            )

            warning_message_pdb_entries_multiple_files.add(pdb_entry_to_entity[key])

    if pdb_entries_in_multiple_files:
        pdbids = []
        entities = []
        chains = []
        cluster_files = []

        for label, cluster_lst in pdb_entries_in_multiple_files.items():
            for cluster in cluster_lst:
                pdbids.append(label[:12])
                chains.append(label[13:])
                entities.append(pdb_entry_to_entity[label])
                cluster_files.append(cluster)

        df = pd.DataFrame({"PDB ID":pdbids,
                      "ENTITY":entities,
                      "CHAIN":chains,
                      "CLUSTER FILE PATH":cluster_files})

        df.to_csv(f"{clustering_directory}/duplicate_entities_warning.csv", index=False)

        print(f"WARNING: There were some entities that contained multiple chains that do not correspond to the same entity. "
    f"We only assign one label type per entity, so one of the chains may be mislabeled in the input file. You can find the affected structures in the following table {clustering_directory}/duplicate_entities_warning.csv."
              f"The table lists the PDB ID, entity ID, chain ID, and corresponding cluster file path for each affected entry.")

    if entities_in_multiple_files:
        entities = []
        cluster_files = []

        for entity, cluster_lst in entities_in_multiple_files.items():
            for cluster in cluster_lst:
                cluster_files.append(cluster)
                entities.append(entity)

        df = pd.DataFrame({"ENTITY": entities,
                      "CLUSTER FILE PATH": cluster_files})

        df.to_csv(f"{clustering_directory}/duplicate_entries_warning.csv", index=False)

        print(
            f"WARNING: There are chains of some entities that may be assigned different chain names in each cluster. "
            "This is common when the dataset contains structures from multiple organisms. To keep one chain name for "
            f"this entity across all clusters, edit the entity list. You can find the affected clusters in the table "
            f"{clustering_directory}/duplicate_entries_warning.csv. "
            "The table lists each affected entity and the cluster file paths in which it appears."
        )


    return cluster_pdbchids, outlier_case_entries


def generate_chain_assignments(current_list, master_molID_class_list, cluster_pdbchids, outlier_case_entries, project_dir):
    """
    Generates consistent chain assignments for polymer and non-polymer entities
    using sequence clustering results and structural entity information.

    Identifies entities and PDB chain IDs appearing across multiple cluster files,
    detects potential conflicts, and assigns chain labels while ensuring that
    identical entities retain the same chain name across different PDB structures.
    Entities not represented in the clustering results are assigned chain labels
    through a fallback mechanism.

    The final entity-to-chain mapping is alphabetically sorted and written to
    an input text file.

    Parameters
    ----------
    current_list (list): List of molID conversion objects containing entity
                         identifiers, occurrence counts, and completeness information.
    master_molID_class_list (list): List of molecular classification objects
                                   containing structural entity information,
                                   including entity names and entity types.
    cluster_pdbchids (dict): Dictionary mapping cluster representatives to their member PDB chain IDs.
    outlier_case_entries (dict): Dictionary mapping outlier PDB entries/entities to their associated labels.
    project_dir (str): Path to the project directory containing the clustering
                       storage directory and where the final input file will
                       be written.

    Returns
    -------
    None
    """

    def get_label(n: int) -> str:
        """
        Generates a unique chain label based on a numerical index.

        Labels are generated sequentially using uppercase letters (A-Z), followed
        by combinations of uppercase letters (AA-ZZ). For indices beyond 701,
        labels are generated using a combination of uppercase letters and numbers.

        Parameters
        ----------
        n (int): Non-negative numerical index used to generate the chain label.
                 Each unique index corresponds to a different chain label.

        Returns
        -------
        str: Generated chain label corresponding to the supplied numerical index.

        Raises
        ------
        ValueError: If the supplied index is negative.
        """
        if n < 0:
            raise ValueError("Input must be a non-negative integer.")

        if n < 702:
            result = []
            n += 1
            while n > 0:
                n, remainder = divmod(n - 1, 26)
                result.append(chr(65 + remainder))
            return "".join(reversed(result))

        # Phase 3: Past 702 -> Letter + Number system
        idx = n - 702
        letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        max_per_letter = 9

        letter_idx = idx // max_per_letter
        num_idx = (idx % max_per_letter) + 1

        letter = letters[letter_idx % 26]
        cycle = (letter_idx // 26) + 1
        if cycle > 1:
            return f"{letter}{num_idx}{cycle}"
        return f"{letter}{num_idx}"

    def get_combination(index):
        """
        Retrieves a specific combination based on its zero-based index.
        Index 0 returns an empty string.
        """
        if index < 0:
            raise IndexError("Index out of range")

        if index == 0:
            return ''

        alphabet = list(string.ascii_lowercase[1:]) + list(string.digits)
        base = len(alphabet)

        length = 1

        while index > base ** length:
            index -= base ** length
            length += 1

        result = ""

        for power in range(length - 1, -1, -1):
            divisor = base ** power
            digit = (index - 1) // divisor
            index = (index - 1) % divisor + 1
            result += alphabet[digit]

        return result

    # 1. INITIALIZATION & DATA EXTRACTION
    clustering_directory = f"{project_dir}/clustering_storage/"
    entity_type_structure = {}
    chain_entity_info = {}
    polymers = {}
    nonpolymers = {}
    conflict_type = []

    for molIDConversion in current_list:
        molIDConversion.check_for_completeness()
        if not molIDConversion.complete:
            chain_count = int(molIDConversion.occur)
            entity_name = str(molIDConversion.molID)
            chain_entity_info[entity_name] = chain_count

    for i in master_molID_class_list:
        polymer_entity = set()
        non_polymer_entity = set()
        structure = os.path.basename(i.file_name)

        for entid in i.entity_type_map:
            entity_name = i.entity_name_map[entid]
            entity_type = i.entity_type_map[entid]
            if entity_type == 'polymer':
                polymer_entity.add(entity_name)
                if entity_name in nonpolymers:
                    conflict_type.append(entity_name)
                else:
                    polymers[entity_name] = 'polymer'
            else:
                non_polymer_entity.add(entity_name)
                if entity_name in polymers:
                    conflict_type.append(entity_name)
                else:
                    nonpolymers[entity_name] = 'nonpolymer'

        entity_type_structure[structure] = {
            'polymer': polymer_entity,
            'nonpolymer': non_polymer_entity
        }

    # 2. SEQUENCE & CLUSTER MAPPINGS
    pdbchid = get_attributes_seqres(f"{clustering_directory}/sub_seqres.txt", attribute='PDBID_CHID')
    entity_names = get_attributes_seqres(f"{clustering_directory}/sub_seqres.txt", attribute='ENTITY_LABEL')

    pdb_entry_to_entity = dict(zip(pdbchid, entity_names))

    # 3. CHAIN LABEL ASSIGNMENT (based on clustering)
    n = 0
    chain_labels = {}  # Maps pdb_chid -> assigned chain letter

    # 3a. Assign ONE label per non-outlier cluster
    for rep, mem_list in cluster_pdbchids.items():

        # Skip clusters containing any outlier PDB chain ID
        is_outlier = any(mem in outlier_case_entries for mem in mem_list)

        if not is_outlier:
            cluster_label = get_label(n)

            for mem in mem_list:
                chain_labels[mem] = cluster_label

            n += 1  # Increment once per cluster

    # 3b. Assign labels to unique outlier entities
    # Tracks the chain label assigned to each entity
    entity_chain_labels = {}

    for entry, entities in outlier_case_entries.items():

        for entity in entities:

            # Assign a new label only if this entity has not been encountered
            if entity not in entity_chain_labels:
                entity_chain_labels[entity] = get_label(n)
                n += 1

            # Reuse the entity's existing label for every PDB chain ID
            if entry not in chain_labels:
                chain_labels[entry] = entity_chain_labels[entity]


    #get_combination
    final_chain_set = {}

    # 4. Build the final entity-to-chains mapping safely
    # We group chains by entity name using the cluster labels we generated
    entity_to_assigned_chains = {}

    # First, collect non-outlier entity assignments
    for entry, chain in chain_labels.items():

        if entry not in outlier_case_entries and entry in pdb_entry_to_entity:
            entity_name = pdb_entry_to_entity[entry]

            if entity_name not in conflict_type:
                entity_to_assigned_chains.setdefault(entity_name, chain)

    # Then, collect outlier entity assignments directly
    # This ensures that multiple entities associated with outliers
    # are not lost when constructing the final mapping.
    for entity, chain in entity_chain_labels.items():

        if entity not in conflict_type:
            entity_to_assigned_chains.setdefault(entity, chain)

    for entity_name, base_chain in entity_to_assigned_chains.items():
        count = chain_entity_info.get(entity_name, 1)
        if entity_name in polymers:
            final_chain_set[entity_name] = [
                f"{base_chain}{get_combination(i)}"
                for i in range(count)
            ]
        else:
            final_chain_set[entity_name] = [base_chain] * count

    # 5. Fallback for any unclustered entities (polymers and non-polymers)
    all_entities = {**polymers, **nonpolymers}
    for entity_name in all_entities:
        if entity_name not in final_chain_set and entity_name in chain_entity_info:
            count = chain_entity_info[entity_name]
            label = get_label(n)
            n += 1

            if entity_name in polymers:
                final_chain_set[entity_name] = [
                    f"{label}{get_combination(i)}"
                    for i in range(count)
                ]
            else:
                final_chain_set[entity_name] = [label] * count

    final_chain_set = dict(sorted(final_chain_set.items()))
    
    # 6. FILE WRITING
    output_file_path = f"{project_dir}/input_file.txt"
    with open(output_file_path, "w") as f:
        for key, val_list in final_chain_set.items():
            joined_values = ",".join(val_list)
            f.write(f"{key}:{joined_values}\n")


def run_clustering_pipeline(filelist, project_dir):
    """
    Runs the initial processing pipeline for preparing structural data
    for sequence clustering.

    Creates the clustering storage directory, extracts sequence information
    from the SEQRES file, identifies the PDB structures of interest from the
    supplied file list, and generates a filtered sub-seqres file containing
    only the selected structures.

    The filtered sequence data is then passed to the RNA and protein
    clustering procedure.

    Parameters
    ----------
    filelist (list): List of paths to the structural files (e.g. CIF files)
                     that will be processed. PDB identifiers are extracted
                     from the file names.
    project_dir (str): Path to the main project directory containing the
                       seqres.txt file and where the clustering storage
                       directory will be created.

    Returns
    -------
    None
    """
    clustering_directory = f"{project_dir}/clustering_storage/"
    seqres_file = f"{project_dir}/seqres.txt"

    try:
        os.mkdir(clustering_directory)
        print(f"Directory '{clustering_directory}' created successfully.\n")

    except FileExistsError:
        print(f"Directory '{clustering_directory}' already exists.")

    seqres_dict = entry_dict_seqres(seqres_file)
    
    PDB_ID_list = []
    
    for file in filelist:
        name = file.split('/')[-1] # pdb_00007qom+01.cif
        pdb_id = (name.split('+')[0]).upper()
        PDB_ID_list.append(pdb_id)

    writing_entity_FASTA(clustering_directory, seqres_dict, PDB_ID_list)

    processing_RNA_PROTEIN_cluster(clustering_directory, f"{clustering_directory}/sub_seqres.txt")