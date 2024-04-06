#!/bin/bash

for i in $(ls *.fasta);

do
phmmer --tblout $i.out $i all_Lpgenomes_prot.faa;
done
