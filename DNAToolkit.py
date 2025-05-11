#DNA Tool-kit file
nucleotide = ['A', 'T', 'C', 'G']

#Checking the sequence is a valid DNA sequence
def validseq(sequence):
    temp_seq = sequence.upper()
    for nuc in temp_seq:
        if nuc not in nucleotide:
            return False
    
    return temp_seq


#Count nucleotides frequence
def count_nuc_frequence(sequence):
    sequence = validseq(sequence)
    if not sequence:
        return "Invalid DNA sequence"
    
    count = {}
    for nuc in nucleotide:
        count[nuc] = sequence.count(nuc)
    
    return count

