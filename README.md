# Lactobacillus Colonization Island Evolution Experiment

Scripts written by Richard Wolff to 
* align reads to WF parent using MIDAS
* perform evolutionary analyses

These scripts were run on UCLA's Hoffman2 computing cluster, using Hoffman2's MIDAS-1.3.2 docker image. 

Documentation:

```
python create_acc_file.py
python convert_gff_to_genes.py

qsub qsub_build_midas_db
qsub qsub_snps_genes
qsub qsub_merge_snps_genes
```

