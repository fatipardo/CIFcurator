#
import os
import shutil
import re
import gzip
import uuid
import numpy as np
from urllib.request import urlopen
from urllib.error import HTTPError, ContentTooShortError
from contextlib import closing, suppress
#
ASSEMBLY_URL='https://files.wwpdb.org/pub/pdb/data/assemblies/mmCIF/all/'
#
def download_pdb_from_metadata(metadata, projdir=None):
    """
    Downloads all the biological assemblies of the PDB entries found in the metadata, in this case
    the description lines of the fasta files, and saves them in the specified project directory.
    Each assembly is saved as its own file: PDBID+NN.cif, where NN is the two-digit assembly number

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
    Downloads a specific PDB file using its ID and saves it in the specified directory.

    Parameters:
    -----------
    pdbid : str
        The PDB ID of the file that will be downloaded.
    pdbformat : str, optional
        The format of the PDB file to be downloaded. Default is '.cif'.
     download_dir : str, optional
        The directory where the downloaded file will be saved. If None, a message will display.

    Returns:
    --------
    None
    """
    download_url='https://files.rcsb.org/download/'
    if download_dir is None:
        print("Please provide a directory where to store downloaded files...")
    else:
        target = download_dir+'/'+pdbid+pdbformat
        source = download_url+pdbid.upper()+pdbformat
        download_from_url(source, target)

def list_assemblies_from_id(pdbid):
    """
    Lists the biological assembly files available in the wwPDB archive for a specific PDB ID.

    The 'all' directory listing is truncated by the server, so the list is read from the
    'divided' directory, which groups the files by the middle two characters of the PDB ID.
    If that listing is truncated too, the assemblies are found by trying assembly1, assembly2, ...

    Parameters:
    -----------
    pdbid : str
        The PDB ID of the entry.

    Returns:
    --------
    filenames : list of str
        The assembly file names (e.g. '1aw1-assembly1.cif.gz'), sorted by assembly number.
    """
    pdbid = pdbid.lower()
    index_url = 'https://files.wwpdb.org/pub/pdb/data/assemblies/mmCIF/divided/'+pdbid[1:3]+'/'
    with closing(urlopen(index_url)) as r:
        listing = r.read().decode('utf-8', errors='replace')
    if 'truncated' not in listing.lower():
        found = re.findall(r'href="(' + pdbid + r'-assembly(\d+)\.cif\.gz)"', listing)
        return [name for name, num in sorted(set(found), key=lambda x: int(x[1]))]
    filenames = []
    while True:
        name = '{0}-assembly{1}.cif.gz'.format(pdbid, len(filenames)+1)
        try:
            with closing(urlopen(ASSEMBLY_URL+name)):
                pass
        except HTTPError as e:
            if e.code == 404:
                return filenames
            raise
        filenames.append(name)

def download_assemblies_from_id(pdbid, download_dir=None):
    """
    Downloads all the biological assemblies of a specific PDB ID from the wwPDB archive,
    decompresses them and saves them in the specified directory as PDBID+NN.cif, where NN is the
    assembly number written with two digits (1aw1-assembly3.cif.gz -> 1aw1+03.cif)

    Parameters:
    -----------
    pdbid : str
        The PDB ID of the entry whose assemblies will be downloaded.
    download_dir : str, optional
        The directory where the downloaded files will be saved. If None, a message will display.

    Returns:
    --------
    None
    """
    if download_dir is None:
        print("Please provide a directory where to store downloaded files...")
    else:
        filenames = list_assemblies_from_id(pdbid)
        if not filenames:
            print('WARNING: no assemblies found for {0}'.format(pdbid))
        for name in filenames:
            number = int(re.search(r'-assembly(\d+)\.cif\.gz$', name).group(1))
            target = download_dir+'/{0}+{1:02d}.cif'.format(pdbid.lower(), number)
            download_from_url(ASSEMBLY_URL+name, target, gunzip=True)

def get_idset_from_metadata(metadata):
    """
    Extracts and returns a set of unique PDB IDs from the provided metadata.

    Parameters:
    -----------
    metadata : list of str
        A list of metadata strings, each containing a PDB ID.

    Returns:
    --------
    idlist : list of str
        A sorted list of unique PDB IDs extracted from the metadata.
    """
    idlist = []
    for elt in metadata:
        idlist.append(elt[1:5])
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
def retrieve_seqfile(seqfile=None):
    """
    Downloads the PDB sequence file from the official RCSB FTP site.

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
    sequrl='https://files.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt'
    if seqfile is None:
        seqfile='seqfile.txt'
    download_from_url(sequrl, seqfile)
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

