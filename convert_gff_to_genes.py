import pandas as pd
import numpy as np

ann_dir = "/u/scratch/r/rwolff/Evolution_Experiment/Annotation"
midas_dir = "/u/scratch/r/rwolff/Evolution_Experiment/midas_db/Lactobacillus"
gff = "LpWF-Nanopore"

df = pd.read_csv(f"{ann_dir}/{gff}.gff",skiprows=1,sep="\t",header=None)

columns = ["gene_id","scaffold_id","start","end","strand","gene_type"]

df = df[[0,2,3,4,6,8]]

df.columns = ["scaffold_id","gene_type","start","end","strand","gene_id"]

df = df[columns]

df["gene_id"] = [g.split("|")[1].split(";")[0] for g in df["gene_id"]]

## Create *.genes and *.mapfile files for MIDAS db build
df.to_csv(f"{midas_dir}/{gff}.genes",index=None,sep="\t")
mapfile = pd.DataFrame(["LpWF-Nanopore","Lactobacillus",1],index=["genome_id","species_id","rep_genome"]).T

## Extract gene descriptions 
df = pd.read_csv(f"{ann_dir}/{gff}.gff",skiprows=1,sep="\t",header=None)
df = df[[0,2,3,4,6,8]]
df.columns = ["scaffold_id","gene_type","start","end","strand","gene_id"]

genes_descriptions = [g.split("|")[1].split(";Name=") for g in df["gene_id"]]
df_desc = pd.DataFrame(pd.Series({g[0]:g[1] for g in genes_descriptions}),columns=["Description"])
df_desc["KEGG Ontology"] = np.nan
for gene in df_desc.index:
    if len(df_desc.loc[gene,"Description"].split(";Ontology_term="))>1:
        df_desc.loc[gene,"KEGG Ontology"] = df_desc.loc[gene,"Description"].split(";Ontology_term=")[1]
        df_desc.loc[gene,"Description"] = df_desc.loc[gene,"Description"].split(";Ontology_term=")[0]
df_desc.to_csv("genes_descriptions.csv",sep="\t")