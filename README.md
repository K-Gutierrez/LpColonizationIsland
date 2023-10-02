# **The genetic basis of *Lactobacillus*-host specificity for the commensal niche in *Drosophila melanogaster* revealed through live imaging of colonization dynamics**

----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

## **1. Paired End Short-Reads Genome assembly**

a) Installing trim-galore and Unicycler

```
conda create --name ShortReadsAssembly
conda install -c bioconda trim-galore
conda install -c "bioconda/label/cf201901" trim-galore

conda install -c bioconda unicycler
conda install -c "bioconda/label/cf201901" unicycle

conda activate ShortReadsAssembly
```

b) Trimming the raw data with Trim-galore

```
mkdir trimgalore-output 

trim_galore -q 28  -o trimgalore-output --path_to_cutadapt My/path/cutadapt --paired --clip_R1 15 --clip_R2 15 --fastqc --dont_gzip *_1.fastq *_2.fastq
```

c) Genome assembly using Unicycler

```
unicycler -1 short_reads_1-trimmed.fastq.gz -2 short_reads_2-trimmed.fastq.gz -o output_dir
```

## **2. HiFi PacBio Genome assembly**

a) Installing SeqTK toolkit, Flye, and Circlator

```
conda create --name HiFiPacBioassembly
conda install -c bioconda seqtk
conda install -c "bioconda/label/cf201901" seqtk

conda install -c bioconda flye
conda install -c "bioconda/label/cf201901" flye

conda install -c bioconda circlator
conda install -c "bioconda/label/cf201901" circulator
```

b) Random sampling 80,000 HiFi PacBio reads using SeqTK toolkit v1.3

```
```



c) HiFi PacBio genome assembly using Flye v2.9.1

```
```
d) circularize HiFi PacBio genome assemblies using Circlator v1.5.5


```
```



```
/data/programs/miniconda3/bin/seqtk sample -s100 /data2/projects/WilliamLudington_QUO1002401/Lactobacillus_plantarum_LpWF_Parental/Lactobacillus_plantarum_LpWF_Parental.hifi_reads.fastq 80000 > /data2/projects/WilliamLudington_QUO1002401/Lactobacillus_plantarum_LpWF_Parental/Lactobacillus_plantarum_LpWF_Parental.subsampled_hifi_reads.fastq
```















