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