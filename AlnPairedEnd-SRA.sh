#Aln Paired-end raw reads 
#!/bin/bash

for i in $(ls My/path/*_1.fastq); do bowtie2 -x dbname -1 $i -2 ${i/_1.fastq/_2.fastq} -S ${i%_1.fastq.gz}.sam --no-unal; 
done
