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

The hybrid assembly of the ancestral strain (3) was used as the reference genome, against which SNVs and CNVs in each replicate/passage were called. 

Documentation:

```
############################################################
Step 1: creating MIDAS database and necessary metadata files
############################################################
python create_acc_file.py
python convert_gff_to_genes.py
qsub qsub_build_midas_db

############################################################################################################
Step 2: Aligning (Illumina short-) reads to reference assembly, comparing SNV and CNV content across samples
############################################################################################################
qsub qsub_snps_genes
qsub qsub_merge_snps_genes

##################################
Step3: create html output document
##################################
jupyter nbconvert --to html --template hidecode snv_gene_trajectories.ipynb

```

*Note*: When performing the inter-sample merge, genes were clustered at the 99% sequence identity threshold (using the ```--cluster_pid 99``` in the ```merge_midas.py``` command of the ```qsub_merge_snps_genes``` script), rather than the default 95%. As in each replicate we are dealing with a simple, single lineage population, this decision was made so that the dynamics of even closely related genes could more easily be discriminated.  
