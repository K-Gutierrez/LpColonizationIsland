# Lactobacillus Colonization Island Evolution Experiment

Scripts written by Richard Wolff to 
* align reads to WF parent using MIDAS
* perform evolutionary analyses

These scripts were run on UCLA's Hoffman2 computing cluster, using Hoffman2's MIDAS-1.3.2 docker image. 

Raw data and annotations furnished by Karina Garcia:

    1) Raw data (Nextseq 2x75) of 12 evolved replicates
    2) Raw data (Nextseq 2X75) for 8 evolved passages
    3) LpWF - Hybrid assembly plus short and long reads (Parental strain that has the colonization island)
    4) LpR3P51 - Hybrid assembly plus short and long reads (Evolved strain that does not have the colonization island)
    5) All the assemblies for the 12 evolved replicates and 8 evolved passages
    
Documentation:

```
python create_acc_file.py
python convert_gff_to_genes.py

qsub qsub_build_midas_db
qsub qsub_snps_genes
qsub qsub_merge_snps_genes
```

