# **The genetic basis of *Lactobacillus*-host specificity for the commensal niche in *Drosophila melanogaster* revealed through live imaging of colonization dynamics**

----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

1. **Paired End Short-Reads Genome assembly**

-Installing trim-galore 

```
conda create --name trimgalore
conda install -c bioconda trim-galore
conda install -c "bioconda/label/cf201901" trim-galore

conda activate trimgalore 
```


- Trimming the raw data with Trim-galore
- 
mkdir trimgalore-output 

trim_galore -q 28  -o trimgalore-output --path_to_cutadapt My/path/cutadapt --paired --clip_R1 15 --clip_R2 15 --fastqc --dont_gzip *_1.fastq *_2.fastq













