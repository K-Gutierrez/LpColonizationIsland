# conda activate biopython
# python ExtractOrderSec.py

# -*- coding: utf-8 -*-
from Bio import SeqIO

# Paths for input and output files
document_path = "/My/path"
fasta_paths = [
    "asp1-gb.fasta",
    "asp2-gb.fasta",
    "asp3-gb.fasta",
    "gftA-gb.fasta",
    "gftB-gb.fasta",
    "secA2-gb.fasta",
    "secY2-gb.fasta"
]
output_path = "/My/Path/salida.fasta"

# Name of the document file
document_file = "Document_A.txt"

# Read IDs from document A
with open(document_path + document_file, "r") as file:
    ids_document_a = file.read().splitlines()

# Create a dictionary to store sequences from the FASTA files
sequences = {}

# Read sequences from the FASTA files
for fasta_file in fasta_paths:
    full_path = document_path + fasta_file
    for record in SeqIO.parse(full_path, "fasta"):
        short_id = record.id[0:7]
        if short_id not in sequences:
            sequences[short_id] = [str(record.seq)]
        else:
            sequences[short_id] += [str(record.seq)]

# Write the sequences to the output file
with open(output_path, "w") as output:
    for id_document_a in ids_document_a:
        output.write(f">{id_document_a}\n")
        for sequences_fasta in sequences.get(id_document_a, []):
            output.write(sequences_fasta)
            output.write("\n")
