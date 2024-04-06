#!/bin/bash

for i in $(ls My/path/*_1.fastq); 
do 
trim_galore -q 25  -o Trimmed --path_to_cutadapt /My/path/trimgalore/env/bin/cutadapt --clip_R1 15 --fastqc --dont_gzip $i;
done
