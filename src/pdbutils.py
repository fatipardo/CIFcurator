#
import os
import shutil
import re
import gzip
import uuid
import json
import numpy as np
from urllib.request import urlopen, Request
from urllib.error import HTTPError, ContentTooShortError
from contextlib import closing, suppress
#
# OCT2026: files come from the wwPDB beta archive, which uses the extended PDB IDs (pdb_00001aw1)
# and will replace the current archive on July 21, 2027 (https://www.wwpdb.org/ftp/pdb-beta-ftp-sites)
ARCHIVE_URL='https://files-beta.wwpdb.org/pub/wwpdb/pdb/'
DOWNLOAD_URL='https://files-beta.wwpdb.org/download/'
# OCT2026: computed structure models (CSMs) are found with the RCSB Search API. RCSB serves their
# coordinates only as BinaryCIF, so the mmCIF file is downloaded from the database that made the model.
# Only AlphaFold DB models are supported for now: ModelArchive models (MA_...) are skipped
RCSB_SEARCH_URL='https://search.rcsb.org/rcsbsearch/v2/query'
RCSB_ENTRY_URL='https://data.rcsb.org/rest/v1/core/entry/'
AFDB_API_URL='https://alphafold.ebi.ac.uk/api/prediction/'
CSM_WARNING = """WARNING: computed structure models (CSMs) from AlphaFold DB contain ONE PROTEIN CHAIN ONLY.
They have no nucleic acids, ligands, waters or other copies of the protein. For example, a Cas9 model has
no guide RNA or target DNA, and a triosephosphate isomerase model is a monomer, while the experimental
structures are dimers. Models cover the full UniProt sequence (signal peptides included), and the
B-factor column holds the model confidence (pLDDT), not a B-factor."""
#
def extended_pdbid(pdbid):
    """
    Returns a PDB ID in the extended format of the wwPDB beta archive: lower case, 'pdb_' followed by
    eight characters. A legacy four-character ID is padded with zeros ('1AW1' -> 'pdb_00001aw1').

    Parameters:
    -----------
    pdbid : str
        A legacy (e.g. '1aw1') or extended (e.g. 'PDB_00001AW1') PDB ID.

    Returns:
    --------
    pdbid : str
        The extended PDB ID, e.g. 'pdb_00001aw1'.

    Raises:
    -------
    ValueError
        If pdbid is neither a legacy nor an extended PDB ID.
    """
    pdbid = pdbid.strip().lower()
    if re.fullmatch(r'[0-9a-z]{4}', pdbid):
        return 'pdb_0000' + pdbid
    if re.fullmatch(r'pdb_[0-9a-z]{8}', pdbid):
        return pdbid
    raise ValueError('not a PDB ID: {0!r}'.format(pdbid))

def entry_url(pdbid):
    """
    Returns the URL of the directory holding all the files of a PDB entry in the wwPDB beta archive.
    Entries are grouped by the two characters before the last one (pdb_00001aw1 -> entries/aw/pdb_00001aw1/).

    Parameters:
    -----------
    pdbid : str
        A legacy or extended PDB ID.

    Returns:
    --------
    url : str
        The URL of the entry directory, ending with '/'.
    """
    pdbid = extended_pdbid(pdbid)
    return ARCHIVE_URL+'data/entries/'+pdbid[-3:-1]+'/'+pdbid+'/'
#
def download_pdb_from_metadata(metadata, projdir=None):
    """
    Downloads all the biological assemblies of the PDB entries found in the metadata, in this case
    the description lines of the fasta files, and saves them in the specified project directory.
    Each assembly is saved as its own file: PDBID+NN.cif, where PDBID is the extended PDB ID
    (pdb_00001aw1) and NN is the two-digit assembly number

    Parameters:
    -----------
    metadata : list of str
        A list of metadata strings, from which PDB IDs will be extracted.
    projdir : str, optional
        The path to the project directory where the PDB files will be saved. If None, a message will display.

    Returns:
    --------
    None
    """
    if projdir is None:
        print("Please provide a project directory ...")
    else:
        download_dir=projdir+'/raw_bank'
        if not os.path.exists(download_dir):
            os.mkdir(download_dir)
        idset = get_idset_from_metadata(metadata)
        for pdbid in idset:
            download_assemblies_from_id(pdbid, download_dir=download_dir)

def download_pdb_from_id(pdbid, pdbformat='.cif', download_dir=None):
    """
    Downloads a specific PDB file using its ID from the wwPDB beta archive, and saves it, decompressed,
    in the specified directory as PDBID.cif, where PDBID is the extended PDB ID (pdb_00001aw1.cif).

    Parameters:
    -----------
    pdbid : str
        The PDB ID (legacy or extended) of the file that will be downloaded.
    pdbformat : str, optional
        The format of the PDB file to be downloaded ('.cif' or '.pdb'). Default is '.cif'.
     download_dir : str, optional
        The directory where the downloaded file will be saved. If None, a message will display.

    Returns:
    --------
    None
    """
    if download_dir is None:
        print("Please provide a directory where to store downloaded files...")
    else:
        pdbid = extended_pdbid(pdbid)
        target = download_dir+'/'+pdbid+pdbformat
        source = DOWNLOAD_URL+pdbid+pdbformat+'.gz'
        download_from_url(source, target, gunzip=True)

def list_assemblies_from_id(pdbid):
    """
    Lists the biological assembly files available in the wwPDB beta archive for a specific PDB ID,
    read from the 'assemblies' directory of the entry.

    Parameters:
    -----------
    pdbid : str
        The PDB ID (legacy or extended) of the entry.

    Returns:
    --------
    filenames : list of str
        The assembly file names (e.g. 'pdb_00001aw1-assembly1.cif.gz'), sorted by assembly number.
        Empty if the entry has no assemblies directory.
    """
    pdbid = extended_pdbid(pdbid)
    try:
        with closing(urlopen(entry_url(pdbid)+'assemblies/')) as r:
            listing = r.read().decode('utf-8', errors='replace')
    except HTTPError as e:
        if e.code == 404:
            return []
        raise
    found = re.findall(r'href="(' + pdbid + r'-assembly(\d+)\.cif\.gz)"', listing)
    return [name for name, num in sorted(set(found), key=lambda x: int(x[1]))]

def download_assemblies_from_id(pdbid, download_dir=None):
    """
    Downloads all the biological assemblies of a specific PDB ID from the wwPDB beta archive,
    decompresses them and saves them in the specified directory as PDBID+NN.cif, where PDBID is the
    extended PDB ID and NN is the assembly number written with two digits
    (pdb_00001aw1-assembly3.cif.gz -> pdb_00001aw1+03.cif)

    Parameters:
    -----------
    pdbid : str
        The PDB ID (legacy or extended) of the entry whose assemblies will be downloaded.
    download_dir : str, optional
        The directory where the downloaded files will be saved. If None, a message will display.

    Returns:
    --------
    None
    """
    if download_dir is None:
        print("Please provide a directory where to store downloaded files...")
    else:
        pdbid = extended_pdbid(pdbid)
        filenames = list_assemblies_from_id(pdbid)
        if not filenames:
            print('WARNING: no assemblies found for {0}'.format(pdbid))
        for name in filenames:
            number = int(re.search(r'-assembly(\d+)\.cif\.gz$', name).group(1))
            target = download_dir+'/{0}+{1:02d}.cif'.format(pdbid, number)
            download_from_url(entry_url(pdbid)+'assemblies/'+name, target, gunzip=True)

def read_json_from_url(url, query=None):
    """
    OCT2026: Reads a JSON document from a URL. If a query is given, it is sent as JSON with a POST request.

    Parameters:
    -----------
    url : str
        The URL to read.
    query : dict, optional
        The JSON query to send. Default is None (GET request).

    Returns:
    --------
    document : dict or list
        The JSON document. None if the server answers 204 (No Content), as the RCSB Search API
        does when nothing is found.
    """
    if query is None:
        request = url
    else:
        request = Request(url, data=json.dumps(query).encode('utf-8'), headers={'Content-Type': 'application/json'})
    with closing(urlopen(request)) as r:
        if getattr(r, 'status', 200) == 204:
            return None
        return json.loads(r.read().decode('utf-8'))

def search_csm_from_keyword(keyword, max_models=None):
    """
    OCT2026: Finds computed structure models (CSMs) in RCSB PDB with a full text search, the same search
    that the RCSB PDB web site does. It is not the search of retrieve_sequence_from_PDB, so the keyword
    can match a different set of proteins.

    Parameters:
    -----------
    keyword : str
        The text to search for, e.g. 'triosephosphate isomerase'.
    max_models : int, optional
        The largest number of CSMs to return, best matches first. Default is None (all of them).

    Returns:
    --------
    csmids : list of str
        The RCSB IDs of the CSMs, e.g. 'AF_AFP60174F1' (AlphaFold DB) or 'MA_MABAKCEPC' (ModelArchive).
    """
    query = {'query': {'type': 'terminal', 'service': 'full_text', 'parameters': {'value': keyword}},
             'return_type': 'entry',
             'request_options': {'results_content_type': ['computational']}}
    if max_models is None:
        query['request_options']['return_all_hits'] = True
    else:
        query['request_options']['paginate'] = {'start': 0, 'rows': int(max_models)}
    result = read_json_from_url(RCSB_SEARCH_URL, query)
    if result is None:
        return []
    return [hit['identifier'] for hit in result.get('result_set', [])]

def csm_cif_url(csmid):
    """
    OCT2026: Returns the URL of the mmCIF file of a computed structure model (CSM), read from the database
    that made it. RCSB PDB says where each CSM comes from (rcsb_comp_model_provenance); for AlphaFold DB
    models, the AlphaFold DB API gives the URL of the latest model version.

    Parameters:
    -----------
    csmid : str
        The RCSB ID of the CSM, e.g. 'AF_AFP60174F1'.

    Returns:
    --------
    url : str
        The URL of the mmCIF file, e.g. 'https://alphafold.ebi.ac.uk/files/AF-P60174-F1-model_v6.cif'.
        None if the model does not come from AlphaFold DB, or AlphaFold DB does not have it.
    """
    provenance = read_json_from_url(RCSB_ENTRY_URL+csmid).get('rcsb_comp_model_provenance', {})
    if provenance.get('source_db') != 'AlphaFoldDB':
        return None
    entry_id = provenance['entry_id']  # e.g. 'AF-P60174-F1'
    try:
        models = read_json_from_url(AFDB_API_URL+entry_id.split('-')[1])
    except HTTPError as e:
        if e.code == 404:
            return None
        raise
    for model in models:
        if model.get('entryId') == entry_id:
            return model['cifUrl']
    return None

def download_csm_from_keyword(keyword, projdir=None, max_models=None):
    """
    OCT2026: Optional step. Finds computed structure models (CSMs) in RCSB PDB with a full text search
    (search_csm_from_keyword) and downloads their mmCIF files to the raw_bank of the project directory,
    next to the experimental structures. Each model is saved as CSMID+01.cif, where CSMID is the RCSB ID
    (AF_AFP60174F1+01.cif), like the first biological assembly of an experimental structure.
    Only AlphaFold DB models are downloaded; other models are skipped with a message.
    CSMs contain one protein chain only (see CSM_WARNING), which is printed first.

    Parameters:
    -----------
    keyword : str
        The text to search for, e.g. 'triosephosphate isomerase'.
    projdir : str, optional
        The path to the project directory where the files will be saved. If None, a message will display.
    max_models : int, optional
        The largest number of CSMs to download, best matches first. Default is None (all of them).

    Returns:
    --------
    csmids : list of str
        The RCSB IDs of the CSMs that were downloaded.
    """
    if projdir is None:
        print("Please provide a project directory ...")
        return []
    print(CSM_WARNING)
    download_dir=projdir+'/raw_bank'
    if not os.path.exists(download_dir):
        os.mkdir(download_dir)
    csmids = search_csm_from_keyword(keyword, max_models=max_models)
    print('{0} computed structure models were found for {1!r}'.format(len(csmids), keyword))
    downloaded = []
    for csmid in csmids:
        url = csm_cif_url(csmid)
        if url is None:
            print('Skipping {0}: no mmCIF file from AlphaFold DB'.format(csmid))
            continue
        download_from_url(url, download_dir+'/{0}+01.cif'.format(csmid))
        downloaded.append(csmid)
    print('{0} computed structure models were downloaded to {1}'.format(len(downloaded), download_dir))
    return downloaded

def get_idset_from_metadata(metadata):
    """
    Extracts and returns a set of unique PDB IDs from the provided metadata.

    Parameters:
    -----------
    metadata : list of str
        A list of metadata strings (fasta description lines), each starting with PDBID_CHAIN,
        e.g. '>pdb_00001aw1_A mol:protein ...'. Legacy IDs ('>1aw1_A ...') are also accepted.

    Returns:
    --------
    idlist : list of str
        A sorted list of unique extended PDB IDs extracted from the metadata.
    """
    idlist = []
    for elt in metadata:
        idlist.append(extended_pdbid(elt.lstrip('>').split()[0].rsplit('_', 1)[0]))
    return sorted(set(idlist))

#
def retrieve_sequence_from_PDB(keyword, mode='sequence', update=True, seqfile=None):
    """
    Retrieves sequences or metadata from a PDB sequence file based on the keyword match.

    Parameters:
    ---------
    keyword : str
        The keyword to search for in the sequence or metadata.
    mode : str, optional
        Specifies whether to match the keyword in the 'sequence' or 'metadata'. Default is 'sequence'.
    update : bool, optional
        If True, the sequence file will be downloaded or updated before searching. Default is True.
    seqfile : str, optional
        The path to the sequence file. If None, the file will be downloaded if update is True.

    Returns:
    --------
    sequence : numpy.ndarray
        A list of sequences that match the keyword.
    metadata : numpy.ndarray
        A list of metadata associated with the matching sequences (fasta files description line).
    """
    if update:
        with suppress(FileNotFoundError):
            os.remove(seqfile) # remove existing seqfile if any
        seqfile = retrieve_seqfile(seqfile=seqfile)
    metadata = []
    sequence = []
    with open(seqfile) as f:
        nextline=False
        prevline='#'
        for line in f:
            if nextline:
                sequence.append(line)
                nextline=False
            else:
                hit = re.findall(keyword, line, flags=re.I)
                if hit:
                    if(mode=='sequence'):
                        metadata.append(prevline)
                        sequence.append(line)
                    elif(mode=='metadata'):
                        metadata.append(line)
                        nextline = True
            prevline=line
    return np.atleast_1d(sequence), np.atleast_1d(metadata)
#
def retrieve_sequences_matching_from_PDB(ref_sequences, seqfile):
    """
    OCT2026: Retrieves the sequences of the PDB sequence file that match any of the reference sequences,
    reading the file only once. It finds the same chains as calling
    retrieve_sequence_from_PDB(seq, mode='sequence', update=False) for every reference sequence: a chain
    matches when its sequence ends with a reference sequence, ignoring case. That reads the whole file once
    per reference sequence, which takes more than an hour for triosephosphate isomerase (642 reference
    sequences). Each matching chain is listed once, in the order of the file.

    Parameters:
    -----------
    ref_sequences : list of str
        The reference sequences, e.g. those returned by retrieve_sequence_from_PDB(keyword, mode='metadata').
    seqfile : str
        The path to the sequence file.

    Returns:
    --------
    sequence : numpy.ndarray
        The sequences that match a reference sequence.
    metadata : numpy.ndarray
        The metadata of the matching sequences (fasta files description line).
    """
    # Reference sequences grouped by length, so the end of each chain is looked up once per length
    refs_by_length = {}
    for seq in ref_sequences:
        seq = seq.strip().upper()
        if seq:
            refs_by_length.setdefault(len(seq), set()).add(seq)
    metadata = []
    sequence = []
    with open(seqfile) as f:
        prevline='#'
        for line in f:
            if not line.startswith('>'):
                seq = line.rstrip('\n').upper()
                for length, refs in refs_by_length.items():
                    if len(seq) >= length and seq[len(seq)-length:] in refs:
                        metadata.append(prevline)
                        sequence.append(line)
                        break
            prevline=line
    return np.atleast_1d(sequence), np.atleast_1d(metadata)
#
def retrieve_seqfile(seqfile=None):
    """
    Downloads the PDB sequence file (pdb_seqres.fasta.gz) from the wwPDB beta archive and decompresses it.

    Parameters:
    -----------
    seqfile : str, optional
        The path where the sequence file will be saved. If None, the file will be saved as 'seqfile.txt'.

    Returns:
    --------
    seqfile : str
        The path to the downloaded sequence file.
    """
    #sequrl='ftp://ftp.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt' FAPA CHANGED ADDRESS FEB 2025
    #sequrl='https://files.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt' OCT2026: beta archive, extended PDB IDs
    sequrl=ARCHIVE_URL+'derived_data/pdb_seqres.fasta.gz'
    if seqfile is None:
        seqfile='seqfile.txt'
    download_from_url(sequrl, seqfile, gunzip=True)
    return seqfile
#
def download_from_url(source, target, gunzip=False):
    """
    Downloads a file from a given URL and saves it to a specified target location.

    Parameters:
    -----------
    source : str
        The URL of the file to be downloaded.
    target : str
        The path where the downloaded file will be saved.
    gunzip : bool, optional
        If True, the downloaded file is gzip-compressed and is decompressed while saving. Default is False.

    Returns:
    --------
    None

    Raises:
    -------
    ContentTooShortError
        If fewer bytes arrive than the server's Content-Length. The data is written to
        a staging file unique to this call, next to target, and renamed only once complete,
        so target is never left truncated, even when several calls download the same target.
    """
    # Unique name, opened with 'x' so it is never shared; unlike tempfile.mkstemp this keeps
    # the usual umask permissions, which the downloaded file inherits.
    partial = '{0}.{1}.part'.format(target, uuid.uuid4().hex)
    # Content-Length counts the compressed bytes, so a gzip file is checked before it is decompressed.
    unzipped = partial + '.gunzip'
    try:
        with closing(urlopen(source)) as r:
            expected = r.info().get('Content-Length')
            with open(partial, 'xb') as f:
                shutil.copyfileobj(r,f)
                received = f.tell()
        if expected is not None and received < int(expected):
            raise ContentTooShortError('retrieval incomplete: got only {0} out of {1} bytes from {2}'.format(
                received, expected, source), None)
        if gunzip:
            with gzip.open(partial, 'rb') as g, open(unzipped, 'xb') as f:
                shutil.copyfileobj(g,f)
            os.replace(unzipped, target)
        else:
            os.replace(partial, target)
    finally:
        for staging in (partial, unzipped):
            with suppress(FileNotFoundError):
                os.remove(staging)
    print('wrote {0} from {1}'.format(target, source))

