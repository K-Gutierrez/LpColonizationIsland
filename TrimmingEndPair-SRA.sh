#For Paired-end
#!/bin/bash

for i in $(ls My/path/*_1.fastq); 
do 
trim_galore -q 25  -o Trimmed --path_to_cutadapt /My/path/trimgalore/env/bin/cutadapt --paired --clip_R1 15 --clip_R2 15 --fastqc --dont_gzip $i ${i/_1.fastq/_2.fastq};
done
