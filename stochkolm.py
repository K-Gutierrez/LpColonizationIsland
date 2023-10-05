#usage 
# conda activate biopython
#python stochkolm.py 

from Bio import SeqIO

records = SeqIO.parse("salida-aln-gb-final.fasta", "fasta")
count = SeqIO.write(records, "aSecSystem", "stockholm")
print("Converted %i records" % count)
