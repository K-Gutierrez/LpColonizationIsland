#Align Single end raw reads (PacBio HiFi)
#!/bin/bash

for i in $(ls My/path/*_1.fastq); 
minimap2 -t 16 --MD --eqx -ayYLx asm20 ColonizationIsland-LpWF-Masked.fasta $i  all_evolved_barcode01.fastq > ${i%_1.fastq}.sam;
done
