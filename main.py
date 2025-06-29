from DNAToolkit import *

sequence = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"

print(f"[1] Valid DNA sequence: {validseq(sequence)}\n")

print(f"[2] number of nucleotide: {count_nuc_frequence(sequence)}\n")

print(f"[3] RNA sequence: {transcription(sequence)}\n")

print(f"[4] DNA string + Compliment + Reverse complement:\n    5' {sequence} 3'")
print(f"       {''.join(['|' for c in range(len(sequence))])}")
print(f"    3' {complement(sequence)} 5' [complement]")
print(f"    5' {reverse_complement(sequence)} 3' [reverse complement]\n")

print(f"[5] mRNA: {transcription(sequence)}\n")

print(f"[6] GC content: {gc_content(sequence)}%\n")
print(f"[7] GC content in subsection k =5: \n    {gc_content_subsec(sequence, k=5)}\n")

seq1 = "GAGCCTACTAACGGGAT"
seq2 = "CATCGTAATGACGGCCT"
print(f"[8] Hamming distance  between {seq1} & {seq2} is: {hamming_distance(seq1, seq2)}\n")

mRNA_sequence = transcription(sequence)
print(f"[9] Protien sequence: {translation(mRNA_sequence)}\n")

print(f"[10] codon usage frequency: {codon_usage(mRNA_sequence)}\n")

amino_acid = "MVLLV"
print(f"[11] Amino acid to mRNA: {amino_acid_to_mRNA(amino_acid)}\n")
 
mRNA = "AUGGCCAUGCUUGAUAG"
intron_positions = [3, 10]  # Insert introns after these indices (0-based)
intron_patterns = [
    ["GTACGTACGT", "TGCATGCATG"],  # Pattern 1
    ["AAAAA", "TTTTT"],             # Pattern 2
    ["GATCGATC", "CCTAGCTA"],       # Pattern 3
    ["CACACACA", "TGTGTGTG"],       # Pattern 4
    ["ATATATAT", "GCGCGCGC"]        # Pattern 5
]

print(f"[12] mRNA to hnRNA: {generate_hnRNA_sequences(mRNA_sequence, intron_positions, intron_patterns, output_file='hnRNA_sequences.fasta')}\n")