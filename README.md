# **The genetic basis of *Lactobacillus*-host specificity for the commensal niche in *Drosophila melanogaster* revealed through live imaging of colonization dynamics**

----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

## **1. Paired End Short-Reads Genome assembly**

a) Installing trim-galore, Unicycler, and Bbmap

```
conda create --name ShortReadsAssembly
conda install -c bioconda trim-galore
conda install -c "bioconda/label/cf201901" trim-galore

conda install -c bioconda unicycler
conda install -c "bioconda/label/cf201901" unicycle

conda install -c bioconda quast
conda install -c "bioconda/label/cf201901" quast

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

d) Genome assembly metrics

```
quast genomeassembly.fasta -o out_name
```



## **2. HiFi PacBio Genome assembly**

a) Installing SeqTK toolkit, Flye, Circlator, Ragoo, and Quast

```
conda create --name HiFiPacBioassembly
conda install -c bioconda seqtk
conda install -c "bioconda/label/cf201901" seqtk

conda install -c bioconda flye
conda install -c "bioconda/label/cf201901" flye

conda install -c bioconda circlator
conda install -c "bioconda/label/cf201901" circulator

conda install -c bioconda quast
conda install -c "bioconda/label/cf201901" quast

conda install -c imperial-college-research-computing ragoo

conda activate HiFiPacBioassembly
```

b) Random sampling 80,000 HiFi PacBio reads using SeqTK toolkit v1.3

```
seqtk sample -s100 HiFiPacBioReads.fastq 80000 > My/path/HiFiPacBioReads.subsampled.fastq
```

c) HiFi PacBio genome assembly using Flye v2.9.1

```
flye -t 20 --out-dir HiFiPacBioAssembly --pacbio-hifi HiFiPacBioReads.subsampled.fastq --genome-size 3.5m

# The genome size of Lactiplantibacillus plantarum is ~3.5 Mbp
```

d) circularize HiFi PacBio genome assemblies using Circlator v1.5.5

```
circlator all HiFiPacBioassembly.fasta HiFiPacBioReads.subsampled.fastq output_directory
```

e) Genome assembly metrics

```
quast HiFiPacBiogenomeassembly.fasta -o out_name
```

f) Ordering contigs using the LpWF-HiFi genome assembly as a reference

```
ragoo.py HiFiPacBio-circlator.fasta LpWF-HiFi.fasta

# The genome assemblies from the R3P51, D11, and F9A were ordered using LpWF-HiFi as a reference
```



## **3. Comparative genomics of L. plantarum genome assemblies and their functional annotation based on RAST subsystems.**

a) 

```

```

b) 

```

```

c) 

```

```

d) 

```

```

e) 

```

```

f) 

```

```

g) 

```

```



## **4. Indels distribution over time in the evolved replicates and passages.**

a) Installing Bowtie2, Samtools, Sambamba, Freebayes, and SnpEff

```
conda create --name Indels

conda install -c bioconda bowtie2
conda install -c "bioconda/label/broken" bowtie2
conda install -c "bioconda/label/cf201901" bowtie2

conda install -c bioconda samtools
conda install -c "bioconda/label/cf201901" samtools

conda install -c bioconda freebayes
conda install -c "bioconda/label/broken" freebayes
conda install -c "bioconda/label/cf201901" freebayes

conda activate Indels

To install Sambamba, visit: https://lomereiter.github.io/sambamba/docs/sambamba-markdup.html
To install SnpEff, visit http://snpeff.sourceforge.net/SnpEff_manual.html#run
```

b) Aligning the short reads to a reference genome

```
# Create a database

bowtie2-build -f ReferenceGenomeAssembly.fasta dbname

# Align the Paired Short-Reads

bowtie2 -x dbname -1 short_reads_1-trimmed.fastq -2 short_reads_2-trimmed.fastq -S ShortReadsAln.sam --no-unal
```

c) Convert SAM to BAM for sorting

```
samtools view -S -b ShortReadsAln.sam > ShortReadsAln.bam
```

d) Sort BAM for SNP calling

```
samtools sort ShortReadsAln.bam ShortReadsAln-sorted.bam
```

e) Mark duplicates with sambamba

```
./sambamba-0.7.1-linux-static markdup -r -p ShortReadsAln-sorted.bam ShortReadsAln-sorted_filter.bam
```

f) Freebayes to find Indels

```
./freebayes-v1.3.1 -f ReferenceGenomeAssembly.fasta --gvcf --use-best-n-alleles 4 -p 1 ShortReadsAln-sorted_filter.bam | vcffilter -f "QUAL > 20" > Indels.vcf
```

g) Split multiallelic variants 

```
bcftools norm --multiallelics -both Indels.vcf > Indels-split.vcf
```

h) Filter only SNPs

```
vcffilter -f "TYPE = snp" Indels-split.vcf > Indels-split-SNPs.vcf
```

i) Annotate the genetic variants using SnpEff

```

```

g) 

```









## **5. Indels distribution over time in the evolved replicates and passages.**

Comparative genomics of Lp WF and Lp R3P51, and the functional annotation of the colonization island.

Mapping Illumina reads on the colonization island
In silico detection of circular and linear contigs
L. plantarum phylogenetic tree.

Genome mining for colonization island in L. plantarum genomes

Transposable elements annotation and classification.
SRRPs and TEs similarity network.
Prediction of the recombination sites.
Colonization island in other bacteria.













