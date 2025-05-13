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

# complement sequence
def complement(sequence):
    complement_seq = ''
    for nuc in sequence:
        complement_seq += DNA_reverse_Complement[nuc]
    
    return complement_seq

