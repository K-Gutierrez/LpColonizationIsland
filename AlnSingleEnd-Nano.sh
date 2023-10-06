#Align Single end raw reads (Nanopore)
#!/bin/bash


for i in $(ls My/path/*_1.fastq); 
minimap2 -t 20 -ax map-ont ColonizationIsland-LpWF-Masked.fasta $i  all_evolved_barcode01.fastq > ${i%_1.fastq}.sam;
done
