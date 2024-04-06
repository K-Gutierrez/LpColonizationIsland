import pandas as pd
import config
import numpy as np
import os

### PacBio settings
ann_dir = f"{config.base_input_dir}/PacBio-Assemblies"

midas_dir = config.midas_dir

gff = "LpParental-HiFi"
    
df = pd.read_csv(f"{ann_dir}/{gff}.gff",skiprows=1,sep="\t",header=None)

columns = ["gene_id","scaffold_id","start","end","strand","gene_type"]

df = df[[0,2,3,4,6,8]]

df.columns = ["scaffold_id","gene_type","start","end","strand","gene_id"]

df = df[columns]

df["gene_id"] = [g.split("|")[1].split(";")[0] for g in df["gene_id"]]

gffmod = gff.split("-")[0]

## Create *.genes and *.mapfile files for MIDAS db build
os.makedirs(f"{midas_dir}/Lactobacillus/{gffmod}",exist_ok=True)
df.to_csv(f"{midas_dir}/Lactobacillus/{gffmod}/{gffmod}.genes",index=None,sep="\t")


## Create mapfile for MIDAS 
M = np.array([["LpParental","Lactobacillus",1]])
mapfile = pd.DataFrame(M.T,index=["genome_id","species_id","rep_genome"]).T
mapfile.to_csv(f"{midas_dir}/Lactobacillus.mapfile",index=None,sep="\t")
