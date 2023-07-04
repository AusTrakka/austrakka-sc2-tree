# Convert the fasta files to tabular format
seqkit fx2tab test.fasta > file1.tab
seqkit fx2tab test2.fasta > file2.tab

# Sort the tabular files
sort file1.tab > file1_sorted.tab
sort file2.tab > file2_sorted.tab

# Compare the sorted files
comm -3 file1_sorted.tab file2_sorted.tab > differences.tab

cut -f1 differences.tab > changed_ids.txt
