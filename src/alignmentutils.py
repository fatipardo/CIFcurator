from __future__ import print_function
from __future__ import division
import os
import subprocess
import tempfile
from functools import lru_cache
from Bio import SeqIO
from Bio import Align
from Bio.Align import substitution_matrices
import numpy as np

# AA Map from 3 letter amino acid id to 1 letter id
# it also includes nucleic acids, including post-tranlational modifications,
# which are mapped to ACTUG.
def ResnConvert(resn):
    """
    Converts the 3 letter amino acid id into a singular letter ID.

    Parameters:
    -----------
    resn : str
        The 3 letter amino acid id

    Returns:
    --------
    ans : str
        The one-letter ID representing the 3 letter amino acid ID. If there's none, "X" is
        returned.

    """
    AA = {}
    AA["UNK"] = "X"
    AA["ALA"] = "A"
    AA["ARG"] = "R"
    AA["ASN"] = "N"
    AA["ASP"] = "D"
    AA["CYS"] = "C"
    AA["GLN"] = "Q"
    AA["GLU"] = "E"
    AA["GLY"] = "G"
    AA["HIS"] = "H"
    AA["ILE"] = "I"
    AA["LEU"] = "L"
    AA["LYS"] = "K"
    AA["MET"] = "M"
    AA["PHE"] = "F"
    AA["PRO"] = "P"
    AA["SER"] = "S"
    AA["THR"] = "T"
    AA["TRP"] = "W"
    AA["TYR"] = "Y"
    AA["VAL"] = "V"
    AA["A"]   = "A"
    AA["C"]   = "C"
    AA["U"]   = "U"
    AA["G"]   = "G"
    AA["2MA"] = "A"
    AA["3AU"] = "U"
    AA["4AC"] = "C"
    AA["4OC"] = "C"
    AA["4SU"] = "U"
    AA["5MC"] = "C"
    AA["5MU"] = "U"
    AA["6IA"] = "A"
    AA["6MZ"] = "U"
    AA["7MG"] = "G"
    AA["8AN"] = "A"
    AA["CM0"] = "C"
    AA["G7M"] = "G"
    AA["H2U"] = "U"
    AA["MIA"] = "A"
    AA["OMC"] = "C"
    AA["OMG"] = "C"
    AA["PSU"] = "U"
    AA["QUO"] = "G"
    AA["T6A"] = "A"
    AA["U8U"] = "U"
    AA["YG"] = "G"
    # NOV2026: DNA (before, every DNA residue was X, so all DNA chains of the same length looked the same).
    # Inosine and unknown nucleotides are N, so the chain is still recognized as a nucleic acid (ACGTUN).
    AA["DA"]  = "A"
    AA["DC"]  = "C"
    AA["DG"]  = "G"
    AA["DT"]  = "T"
    AA["DU"]  = "U"
    AA["DI"]  = "N"
    AA["DN"]  = "N"

    if resn not in AA:
        ans = "X"
    else:
        ans = AA[resn]
    return ans
# END AA Map from 3 letter amino acid id to 1 letter id

# NOV2026: pairwise alignments (Step 4) are done by Biopython, in this process, instead of starting MUSCLE
# for every pair of chains: about 1 ms instead of about 1 s per alignment, with the same identity scores
# (within 0.005 on TIM chains). Global alignment with BLOSUM62 and the gap scores of EMBOSS needle;
# nucleic acid chains (only ACGTUN, see ResnConvert) are scored with matches and mismatches.
PROTEIN_MATRIX = substitution_matrices.load('BLOSUM62')
NUCLEIC_LETTERS = set('ACGTUN')

def make_pairwise_aligner(nucleic=False):
    """
    NOV2026: Creates the Biopython aligner used for pairwise alignments.

    Parameters:
    -----------
    nucleic : bool, optional
        If True, the aligner is for nucleic acid sequences (match 5, mismatch -4, as EMBOSS needle with
        EDNAFULL). If False (default), for proteins (BLOSUM62).

    Returns:
    --------
    aligner : Bio.Align.PairwiseAligner
        A global aligner with gap open score -10 and gap extension score -0.5.
    """
    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    if nucleic:
        aligner.match_score = 5
        aligner.mismatch_score = -4
    else:
        aligner.substitution_matrix = PROTEIN_MATRIX
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    return aligner

PROTEIN_ALIGNER = make_pairwise_aligner()
NUCLEIC_ALIGNER = make_pairwise_aligner(nucleic=True)

def AlignPair(seq1, seq2):
    """
    NOV2026: Aligns two sequences with Biopython (see make_pairwise_aligner) and returns them with gaps ('-').
    Letters that BLOSUM62 does not know are aligned as X, and given back unchanged.

    Parameters:
    -----------
    seq1 : str
        The first sequence (one-letter codes, as from ResnConvert).
    seq2 : str
        The second sequence.

    Returns:
    --------
    aligned_seq : list of str
        The two sequences, aligned: both have the same length.
    """
    seq1, seq2 = seq1.upper(), seq2.upper()
    if not seq1 or not seq2:
        length = max(len(seq1), len(seq2))
        return [seq1.ljust(length, '-'), seq2.ljust(length, '-')]
    if set(seq1) <= NUCLEIC_LETTERS and set(seq2) <= NUCLEIC_LETTERS:
        aligner, known = NUCLEIC_ALIGNER, NUCLEIC_LETTERS
    else:
        aligner, known = PROTEIN_ALIGNER, set(PROTEIN_MATRIX.alphabet)
    masked = [''.join(c if c in known else 'X' for c in seq) for seq in (seq1, seq2)]
    alignment = aligner.align(masked[0], masked[1])[0]
    aligned_seq = []
    for seq, row in zip((seq1, seq2), (alignment[0], alignment[1])):
        letters = iter(seq)
        aligned_seq.append(''.join('-' if c == '-' else next(letters) for c in row))
    return aligned_seq

@lru_cache(maxsize=None)
def PairwiseIdentity(seq1, seq2):
    """
    NOV2026: Aligns two sequences (AlignPair) and returns their identity score (ScoreSequenceAlignment).
    Scores are remembered, so the many identical chains of an ensemble are only aligned once.

    Parameters:
    -----------
    seq1 : str
        The reference sequence, e.g. a standard sequence.
    seq2 : str
        The sequence being compared, e.g. a chain of a structure.

    Returns:
    --------
    score : float
        The fraction of the alignment columns where both sequences have the same letter.
    """
    aligned_seq = AlignPair(seq1, seq2)
    return ScoreSequenceAlignment(aligned_seq[0], aligned_seq[1])

def run_muscle(args):
    """
    NOV2026: Runs MUSCLE and waits until it is done (before: os.popen and checks every second, or every
    10 seconds, for the output file). Raises an error if MUSCLE fails, instead of waiting forever.

    Parameters:
    -----------
    args : list of str
        The arguments of muscle, e.g. ['-align', 'in.fa', '-output', 'out.afa'].

    Returns:
    --------
    None
    """
    try:
        subprocess.run(['muscle'] + list(args), check=True)
    except FileNotFoundError:
        raise RuntimeError('muscle was not found: install it (e.g. conda install -c bioconda muscle) '
                           'and make sure it is on the PATH')

# NOV2026: muscle -super5 stops with "usorter.cpp(88) assert failed: L >= m_WordLength" when a sequence
# is shorter than 8 residues (e.g. a DNA strand with only a few nucleotides modelled); -align does not
SUPER5_MIN_LENGTH = 8

def run_muscle_msa(sequence_vec, seq_ids, fafile, output):
    """
    NOV2026: Writes the multiple sequence alignment of sequence_vec to output, in the format of MUSCLE
    ('> Seq <id>' headers). Uses muscle -super5, or muscle -align when a sequence is shorter than
    SUPER5_MIN_LENGTH. MUSCLE crashes on empty sequences and on a single sequence, so empty sequences
    are left out of MUSCLE and written as all gaps, and with fewer than two non-empty sequences MUSCLE
    is not run at all.

    Parameters:
    -----------
    sequence_vec : list of str
        The sequences to align, as written in fafile.
    seq_ids : list of str
        The identifier of each sequence (header '> Seq <id>').
    fafile : str
        The FASTA file with all the sequences (input of MUSCLE when no sequence is empty).
    output : str
        The aligned FASTA file to write.

    Returns:
    --------
    None
    """
    nonempty = [(seq_id, seq) for seq_id, seq in zip(seq_ids, sequence_vec) if len(seq) > 0]
    empty_ids = [seq_id for seq_id, seq in zip(seq_ids, sequence_vec) if len(seq) == 0]

    if len(nonempty) >= 2:
        mode = '-super5' if min(len(seq) for _, seq in nonempty) >= SUPER5_MIN_LENGTH else '-align'
        if mode == '-align':
            print("Sequences shorter than " + str(SUPER5_MIN_LENGTH) + " residues in " + fafile
                  + ", so I will use muscle -align instead of -super5")
        if empty_ids:
            muscle_input = output + ".in.fa"
            with open(muscle_input, 'w') as f:
                for seq_id, seq in nonempty:
                    f.write("> Seq " + str(seq_id) + "\n" + seq + "\n")
        else:
            muscle_input = fafile
        try:
            run_muscle([mode, muscle_input, '-output', output])
        finally:
            if muscle_input != fafile and os.path.exists(muscle_input):
                os.remove(muscle_input)
        with open(output) as f:
            aligned_length = sum(len(line.strip()) for line in f.read().split('>')[1].split('\n')[1:])
    else:
        with open(output, 'w') as f:
            for seq_id, seq in nonempty:
                f.write("> Seq " + str(seq_id) + "\n" + seq + "\n")
        aligned_length = len(nonempty[0][1]) if nonempty else 0

    if empty_ids:
        print("Empty sequences in " + fafile + ", written as gaps: " + ", ".join(str(i) for i in empty_ids))
        with open(output, 'a') as f:
            for seq_id in empty_ids:
                f.write("> Seq " + str(seq_id) + "\n" + "-" * aligned_length + "\n")

def AlignSequences(sequence_vec):
    """
    Takes a list of sequence strings and performs a MUSCLE alignment,
    outputting a vector of aligned sequence strings.
    NOV2026: two sequences are aligned with Biopython (AlignPair), more than two with MUSCLE,
    in a temporary directory (before: Seq.fa and Seq.afa in the current directory).

    Parameters:
    -----------
    sequence_vec : list of str
        list of sequences to be aligned

    Returns:
    --------
    aligned_seq : list of str
        A list containing the aligned sequences
    """
    if len(sequence_vec) == 2:
        return AlignPair(sequence_vec[0], sequence_vec[1])

    with tempfile.TemporaryDirectory(prefix='cifcurator_') as tmpdir:
        fafile = os.path.join(tmpdir, 'Seq.fa')
        afafile = os.path.join(tmpdir, 'Seq.afa')
        with open(fafile, 'w') as newfafile:
            for i, seq in enumerate(sequence_vec):
                newfafile.write("> Seq" + str(i) + "\n") # NOV2026: numbered, MUSCLE reorders the sequences
                newfafile.write(seq + "\n")

        run_muscle(['-align', fafile, '-output', afafile])

        aligned_seq_map = {}
        with open(afafile) as seqfile:
            key = None
            for line in seqfile:
                if (line[0] == ">"):
                    key = int(line[1:].strip()[3:])
                    aligned_seq_map[key] = ""
                else:
                    aligned_seq_map[key] += line.strip()
    return [aligned_seq_map[i] for i in range(len(sequence_vec))]

def case_safe_name(file_name):
    """
    Return file_name with an underscore before every lowercase letter of its last part, e.g.
    'Ab' -> 'A_b' and 'AB' -> 'AB'. NOV2026: macOS (and Windows) file names ignore case, so the
    alignments of chains AB and Ab were both written to (and read from) the same AB.fasta.
    Chain IDs have only letters and digits, so two chain IDs never give the same name.
    """
    folder, name = os.path.split(file_name)
    return os.path.join(folder, ''.join('_' + c if c.islower() else c for c in name))


def AlignSequenceResNum(sequence_vec, file_name, this_chainsseq_list_ids):
    """
    Takes a list of sequence strings and performs a MUSCLE alignment, outputting
    a vector of aligned sequence strings. This version checks if an alignment has
    already been provided, before running muscle, and in that case just reads the
    existing alignment.

    Parameters:
    -----------
    sequence_vec : list of str
        list containing sequences from FASTA files
    file_name : str
        Name given to FASTA file (the chain ID). An underscore is added before every lowercase
        letter (see case_safe_name), so chain Ab uses A_b.fa and A_b.fasta.
    this_chainsseq_list_ids : list of str
        List of identifiers for each sequence which are also used as headers in the
        FASTA file.

    Returns:
    --------
    aligned_seq_map : dict
        A dictionary where the keys are the sequence identifiers from `this_chainsseq_list_ids`
        and the values are the corresponding aligned sequence strings
    gap_percentages : np.ndarray
        An array where each element represents the percentage of gaps at that position
        across all sequences.
    """
    chain_name = file_name
    file_name = case_safe_name(file_name)
    if os.path.exists(file_name+".fasta") == False:

        with open(file_name+".fa", 'w') as newfafile:
            i = 0
            for seq in sequence_vec:
                newfafile.write("> Seq " + str(this_chainsseq_list_ids[i]) + "\n")
                newfafile.write(seq + "\n")
                i += 1

        # FAPA 26AUG2025: muscle -super5. NOV2026: waits for MUSCLE to finish, and writes to a temporary
        # file renamed at the end, so an interrupted run never leaves a partial alignment to be reused
        partial = file_name+".fasta.part"
        run_muscle_msa(sequence_vec, this_chainsseq_list_ids, file_name+".fa", partial)
        os.replace(partial, file_name+".fasta")

    else:
        print("Alignment already exists, so I will use that one!")

    aligned_seq_map = {}
    aligned_seq = []
    seq = ""
    with open(file_name + ".fasta") as seqfile:
        for line in seqfile:
            if (line[0] == ">"):
                # Very first line
                if (seq == ""):
                    line = line.strip()
                    line = line.split()
                    key = line[2]
                else:
                    aligned_seq_map[key] = seq
                    seq = ""
                    line = line.strip()
                    line = line.split()
                    key = line[2]
            else:
                seq += line.strip()
        aligned_seq_map[key] = seq

    for item in (aligned_seq_map.keys()):
        aligned_seq.append(aligned_seq_map[item])

    print(chain_name)

    sequences = read_fasta_files( file_name + ".fasta")
    gap_percentages = calculate_gap_percentages(sequences)

    print("Gap percentages per position:")
    print(gap_percentages)

    return (aligned_seq_map,gap_percentages)

# The functions below are used to calculate the percentage of gaps per position
def read_fasta_files(fasta_file):
    """
    Reads FASTA files and extracts the sequences into a list.

    Parameters:
    -----------
    fasta_file : str
        Path to FASTA file containing all the sequences

    Returns:
    --------
    sequences : list of str
        list of sequences from a FASTA file
    """
    sequences = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences.append(str(record.seq))
    return sequences

def calculate_gap_percentages(sequences):
    """
    Calculates the percentage of gaps at each position in a list of sequences.

    Parameters:
    -----------
    sequences : list of str
        list of sequences extracted from a FASTA file

    Returns:
    --------
    gap_percentages : np.ndarray
        An array where each element represents the percentage of gaps at that position
        across all sequences.
    """
    sequence_length = len(sequences[0])
    gap_counts = np.zeros(sequence_length)

    for sequence in sequences:
        for i, char in enumerate(sequence):
            if char == '-':
                gap_counts[i] += 1

    total_sequences = len(sequences)
    gap_percentages = (gap_counts / total_sequences) * 100
    return gap_percentages

# FAPA JULY TEST ENDS

def ScoreSequenceAlignment(seq1, seq2):
    """
    Compares the reference sequence to another sequence and counts for similarity based on
    exact matches between corresponding elements from the two sequences.

    Parameters:
    -----------
    seq1 : str
        The reference sequence
    seq2 : str
        The sequence being compared

    Returns:
    --------
    score : float
        The similarity score between the reference sequence and the sequence being compared.
    """
    # Scores based on exact identity. Should  maybe be updated to take longer
    # of sequences so that it can be used with unaligned seq strings too
    score = 0
    for i in range(len(seq1)):
        if (seq1[i] == seq2[i]):
            score += 1
    score = score/len(seq1)
    return score
