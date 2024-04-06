#Align Single end raw reads 
#!/bin/bash
for i in $(ls My/path/*_1.fastq); do bowtie2 -x dbname -U $i -S ${i%_1.fastq}.sam --no-unal; 
done
