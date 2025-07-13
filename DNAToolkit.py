#DNA Tool-kit file

from structures import *

sequence = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"

nucleotide = ['A', 'T', 'C', 'G']
DNA_reverse_Complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}


#Checking the sequence is a valid DNA sequence
def validseq(sequence):
    temp_seq = sequence.upper()
    for nuc in temp_seq:
        if nuc not in nucleotide:
            return False
    
    return temp_seq


#Count nucleotides frequence
def count_nuc_frequence(sequence):
    
    count = {}
    for nuc in nucleotide:
        count[nuc] = sequence.count(nuc)
    
    return count

#transcription
def transcription(sequence):
    
    RNA_sequence = sequence.replace('T', 'U')
    return RNA_sequence

# complement function
def complement(sequence):
    complement_seq = ''
    for nuc in sequence:
        complement_seq += DNA_reverse_Complement[nuc]
    
    return complement_seq

# reverse complement function
def reverse_complement(sequence):
    reverse_seq = sequence[::-1]
    reverse_seq = complement(reverse_seq)

    return reverse_seq

# GC-content calculator function
def gc_content(sequence):
    gc_content = ((sequence.count('G') + sequence.count('C'))/len(sequence)*100)
    return gc_content

# GC-content subset
'''GC content in DNA/RNA sub-sequence length k, k=20 by default'''
def gc_content_subsec(sequence, k=20):
    res = []
    for i in range(0, len(sequence) -k + 1, k):
        subseq = sequence[i:i+k]
        res.append(gc_content(subseq))
    return res

# Hamming distance function 
'''for counting point mutations between two sequences'''
def hamming_distance(seq1, seq2):
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must be of the same length")
    
    distance = 0
    for nuc1, nuc2 in zip(seq1, seq2):
        if nuc1 != nuc2:
            distance += 1
    return distance

#translation function
def translation(sequence):
    codon_table = {
        'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'AGA': 'R', 'AGG': 'R',
        'AAU': 'N', 'AAC': 'N',
        'GAU': 'D', 'GAC': 'D',
        'UGU': 'C', 'UGC': 'C',
        'GAA': 'E', 'GAG': 'E',
        'CAA': 'Q', 'CAG': 'Q',
        'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
        'CAU': 'H', 'CAC': 'H',
        'AUU': 'I', 'AUC': 'I', 'AUA': 'I',
        'UUA': 'L', 'UUG': 'L', 'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
        'AAA': 'K', 'AAG': 'K',
        'AUG': 'M',  # Start codon
        'UUU': 'F', 'UUC': 'F',
        'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S', 'AGU': 'S', 'AGC': 'S',
        'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'UGG': 'W',
        'UAU': 'Y', 'UAC': 'Y',
        'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
        'UAA': '', 'UAG': '', 'UGA': ''  # Stop codons
    }
    
    sequence = sequence.replace('T', 'U')  # Convert DNA to RNA
    protein = ''
    
    for i in range(0, len(sequence) - 2, 3):
        codon = sequence[i:i+3]
        if codon in codon_table:
            protein += codon_table[codon]
            if codon_table[codon] == '':
                break
    
    return protein

# Finding a Motif in a DNA sequence
def find_substring_locations(s: str, t: str) -> list[int]:
    """
    Finds all locations of string t as a substring of string s.

    Args:
        s: The main DNA string.
        t: The substring to find.

    Returns:
        A list of 1-indexed starting positions of t in s.
    """
    locations = []
    len_s = len(s)
    len_t = len(t)

    # Iterate through s up to the point where t can still fit
    for i in range(len_s - len_t + 1):
        # Extract a substring of s with the same length as t
        current_substring = s[i : i + len_t]

        # Compare the extracted substring with t
        if current_substring == t:
            # If they match, record the 1-indexed starting position
            locations.append(i + 1)

    return locations

def codon_usage(sequence):
    """
    Calculate the frequency of each codon in a given DNA sequence.

    Args:
        sequence: A valid DNA sequence.

    Returns:
        A dictionary with codons as keys and their frequencies as values.
    """
    sequence = transcription(sequence)  # Convert to RNA
    codon_freq = {}
    
    for i in range(0, len(sequence) - 2, 3):
        codon = sequence[i:i+3]
        if len(codon) == 3:  # Ensure we have a complete codon
            if codon in codon_freq:
                codon_freq[codon] += 1
            else:
                codon_freq[codon] = 1
    
    return codon_freq

# function to amino acid to mRNA:
def amino_acid_to_mRNA(amino_acid):
    # Example: Convert amino acid sequence to a possible mRNA sequence

    # Minimal codon table (amino acid: [codon1, codon2, ...])
    codon_table = {
    'A': ['GCU', 'GCC', 'GCA', 'GCG'],
    'R': ['CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
    'N': ['AAU', 'AAC'],
    'D': ['GAU', 'GAC'],
    'C': ['UGU', 'UGC'],
    'Q': ['CAA', 'CAG'],
    'E': ['GAA', 'GAG'],
    'G': ['GGU', 'GGC', 'GGA', 'GGG'],
    'H': ['CAU', 'CAC'],
    'I': ['AUU', 'AUC', 'AUA'],
    'L': ['UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG'],
    'K': ['AAA', 'AAG'],
    'M': ['AUG'],
    'F': ['UUU', 'UUC'],
    'P': ['CCU', 'CCC', 'CCA', 'CCG'],
    'S': ['UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC'],
    'T': ['ACU', 'ACC', 'ACA', 'ACG'],
    'W': ['UGG'],
    'Y': ['UAU', 'UAC'],
    'V': ['GUU', 'GUC', 'GUA', 'GUG'],
    '*': ['UAA', 'UAG', 'UGA']  # Stop codons
    }


    mrna = []
    for aa in amino_acid:
        if aa in codon_table:
            # Choose the first codon for each amino acid
            mrna.append(codon_table[aa][0])
        else:
            raise ValueError(f"Unrecognized amino acid: {aa}")
    return ''.join(mrna)

# mRNA to hnRNA with introns
def generate_hnRNA_sequences(mrna, intron_positions, intron_patterns, output_file=None):
    """
    Generate multiple hnRNA sequences by inserting different intron patterns into mRNA.
    
    Parameters:
    - mrna: mRNA sequence (string)
    - intron_positions: List of positions (indices) to insert introns
    - intron_patterns: List of lists, each containing intron sequences for one pattern
    - output_file: (Optional) Filename to save FASTA sequences. If None, no file is saved.
    
    Returns:
    - List of tuples: [(name, sequence), ...]
    """
    def add_introns(mrna, intron_positions, intron_sequences):
        """Helper: Insert introns into mRNA at specified positions."""
        if len(intron_positions) != len(intron_sequences):
            raise ValueError("Number of intron positions must match number of intron sequences.")
        sorted_pairs = sorted(zip(intron_positions, intron_sequences), key=lambda x: x[0], reverse=True)
        hnRNA = mrna
        for pos, intron in sorted_pairs:
            hnRNA = hnRNA[:pos] + intron + hnRNA[pos:]
        return hnRNA

    hnRNA_sequences = []
    for i, pattern in enumerate(intron_patterns):
        hnRNA = add_introns(mrna, intron_positions, pattern)
        hnRNA_sequences.append((f"hnRNA_pattern_{i+1}", hnRNA))

    if output_file:
        with open(output_file, "w") as f:
            for name, seq in hnRNA_sequences:
                f.write(f">{name}\n{seq}\n")

    return hnRNA_sequences

# functuon for parse the PDB file
def parse_pdb_file(pdb_file):
    """
    Parse a PDB file and extract atom information.
    
    Args:
        pdb_file (str): Path to the PDB file.
    
    Returns:
        list: A list of Atom objects containing atom information.
    """
    atoms = []
    with open(pdb_file, 'r') as file:
        for line in file:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                atom = {
                    "atom_serial": int(line[6:11].strip()),
                    "atom_name": line[12:16].strip(),
                    "residue_name": line[17:20].strip(),
                    "chain_id": line[21].strip(),
                    "residue_sequence_number": int(line[22:26].strip()),
                    "x": float(line[30:38].strip()),
                    "y": float(line[38:46].strip()),
                    "z": float(line[46:54].strip()),
                    "element": line[76:78].strip(),
                    "occupancy": float(line[54:60].strip()),
                }
                atoms.append(atom)
    return atoms
