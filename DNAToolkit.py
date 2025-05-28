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
        'AUG': 'M', 'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
        'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L', 'AUU': 'I',
        'AUC': 'I', 'AUA': 'I', 'GUU': 'V', 'GUC': 'V', 'GUA': 'V',
        'GUG': 'V', 'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 
        # Add remaining codons...
    }
    
    sequence = sequence.replace('T', 'U')  # Convert DNA to RNA
    protein = ''
    
    protein_list = []
    
    for i in range(0, len(sequence) - 2, 3):
        codon = sequence[i:i+3]
        if codon in codon_table:
            protein += codon_table[codon]
    
    return protein
