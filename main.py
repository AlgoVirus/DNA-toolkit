from DNAToolkit import *

sequence = "ATCGATCGATCG"

print("[1] Valid DNA sequence:", validseq(sequence))

print("[2] number of nucleotide:", count_nuc_frequence(sequence))

print("[3] RNA sequence:", transcription(sequence))

print(f"[4] DNA string + Reverse complement:\n 5' {sequence} 3'")
print(f"    {''.join(['|' for c in range(len(sequence))])}")
print(f" 3' {reverse_complement(sequence)} 5'")


