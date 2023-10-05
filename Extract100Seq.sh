#!/bin/bash

# Name of the fasta file
fasta_file="all_SelectedGenomes.faa"

# File with the IDs
ids_file="test.txt"

# Read the IDs from the file
readarray -t SecA_Ids < "$ids_file"

# Iterate over the IDs
for id in "${SecA_Ids[@]}"
do
    # Search for the ID in the fasta file and extract 50 sequences before and after
    grep -A 800 -B 800 "$id" "$fasta_file" >> results.txt
    echo "ID: $id" >> results.txt
done
