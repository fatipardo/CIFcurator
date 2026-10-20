from __future__ import print_function
from __future__ import division
import os
import time
from Bio import SeqIO
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

    if resn not in AA:
        ans = "X"
    else:
        ans = AA[resn]
    return ans
# END AA Map from 3 letter amino acid id to 1 letter id

def AlignSequences(sequence_vec, file_name='Seq', this_chainsseq_list_ids=None, gap_percentages=False):
    """
    Takes a list of sequence strings and performs a MUSCLE alignment,
    outputting a vector of aligned sequence strings.

    Parameters:
    -----------
    sequence_vec : list of str
        list of sequences to be aligned

    Returns:
    --------
    aligned_seq : list of str
        A list containing the aligned sequences
    """
    if this_chainsseq_list_ids is None:
        this_chainsseq_list_ids = [str(i) for i in range(len(sequence_vec))]

    fasta_path = file_name + ".fasta"
    fa_path = file_name + ".fa"

    if not os.path.exists(fasta_path):
        with open(fa_path, 'w') as newfafile:
            for i, seq in enumerate(sequence_vec):
                newfafile.write(f"> Seq {this_chainsseq_list_ids[i]}\n")
                newfafile.write(seq + "\n")

        command = f"muscle -super5 {fa_path} -output {fasta_path}"
        process = os.popen(command)

        while not os.path.exists(fasta_path):
            time.sleep(1)
            print("waiting...")

        while os.path.getsize(fasta_path) == 0:
            time.sleep(1)
            print("waiting even more...")

        process.close()
    else:
        print(f"Alignment already exists {fasta_path}, so I will use that one!")

    aligned_seq_map = {}
    with open(fasta_path) as seqfile:
        key = ""
        seq = ""
        for line in seqfile:
            if line.startswith(">"):
                if seq != "":
                    aligned_seq_map[key] = seq
                    seq = ""
                line = line.strip().split()
                key = line[2] if len(line) > 2 else line[0]
            else:
                seq += line.strip()
        if key:
            aligned_seq_map[key] = seq

    aligned_seq = list(aligned_seq_map.values())

    if gap_percentages:
        sequences = read_fasta_files(file_name + ".fasta")
        gap_percentages = calculate_gap_percentages(sequences)

        print("Gap percentages per position:")
        print(gap_percentages)

        return (aligned_seq_map, gap_percentages)

    else: return aligned_seq

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
