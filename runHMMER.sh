#!/bin/bash

for i in $(ls *.fasta);

do
phmmer --tblout $i.out $i all_lpseq.faa;
done
