# Read mapping count
#!/bin/bash 

for i in $(ls *.bam);

do
htseq-count -a 2 -r pos -t CDS -f bam $i Parental_island_nomobileelem.gff > $i.final_count;
done

