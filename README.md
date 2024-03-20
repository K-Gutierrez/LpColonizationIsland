# **A conserved genetic basis for commensal-host specificity through live imaging of colonization dynamics**

----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

## **1. Paired End Short-Reads Genome assembly**

a) Installing trim-galore, Unicycler, and Bbmap

```
conda create --name ShortReadsAssembly
conda activate ShortReadsAssembly

conda install -c bioconda trim-galore
conda install -c "bioconda/label/cf201901" trim-galore

conda install -c bioconda unicycler
conda install -c "bioconda/label/cf201901" unicycle

conda install -c bioconda quast
conda install -c "bioconda/label/cf201901" quast

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
conda activate HiFiPacBioassembly

conda install -c bioconda seqtk
conda install -c "bioconda/label/cf201901" seqtk

conda install -c bioconda flye
conda install -c "bioconda/label/cf201901" flye

conda install -c bioconda circlator
conda install -c "bioconda/label/cf201901" circulator

conda install -c bioconda quast
conda install -c "bioconda/label/cf201901" quast

conda install -c imperial-college-research-computing ragoo

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



## **3. Comparative genomics of *L. plantarum* genome assemblies and their functional annotation based on RAST subsystems.**

```
The comparison was performed using the Function-based comparison tool from RAST
```


## **4. Indels distribution over time in the evolved replicates and passages.**

a) Installing Bowtie2, Samtools, Sambamba, Freebayes, and SnpEff

```
conda create --name Indels
conda activate Indels

conda install -c bioconda bowtie2
conda install -c "bioconda/label/broken" bowtie2
conda install -c "bioconda/label/cf201901" bowtie2

conda install -c bioconda samtools
conda install -c "bioconda/label/cf201901" samtools

conda install -c bioconda freebayes
conda install -c "bioconda/label/broken" freebayes
conda install -c "bioconda/label/cf201901" freebayes

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
# Modify the snpEff.config file, adding the following information

#Lactobacillus plantarum, tig00000001_pilon / Assembly_4Iterations
lp.genome : Lactobacillus_plantarum
lp.chromosomes : tig00000001_pilon / Assembly_4Iterations
lp.tig00000001_pilon.codonTable: Bacterial_and_Plant_Plastid

# Generate the Lp database

java -jar snpEff.jar build -gff3 -v lp

# Run the annotation

java -Xmx4G -jar snpEff.jar lp Indels-split-SNPs.vcf > Indels-split-SNPs-Annotation.vcf
```



## **5. Comparative genomics of Lp WF and Lp R3P51, and the functional annotation of the colonization island.**

a) Blastn using the LpWF as a reference genome

```
makeblastdb -in LpR3P51.fasta  -dbtype nucl

blastall -p blastn -i A.fasta -d B.fasta -m 6 > WF-R3P51Blast.tsv
```

b) Genome Annotation using Interproscan 

```
interproscan.sh -i GenomeAssembly.fasta -f tsv -iprlookup -goterms -pa
```


## **6. Mapping Illumina reads on the colonization island.**

a) Installing Bowtie2, Bedtools, bbmap, and Picard.

```
conda create --name MappingIlluminaReads
conda activate MappingIlluminaReads

conda install -c bioconda bowtie2
conda install -c "bioconda/label/broken" bowtie2
conda install -c "bioconda/label/cf201901" bowtie2

conda install -c bioconda bbmap
conda install -c "bioconda/label/cf201901" bbmap

To install Bedtools, visit: https://bedtools.readthedocs.io/en/latest/content/installation.html
To install Picard, visit: https://broadinstitute.github.io/picard/
```

b) Masking all the TEs within the colonization island from LpWF

```
bedtools maskfasta -fi ColonizationIsland-LpWF.fasta -bed ColonizationIsland-LpWF.bed  -fo ColonizationIsland-LpWF-Masked.fasta
```

c)  Normalizing the Illumina short reads to 100x 

```
/mypath/bbnorm.sh in=short_reads_1-trimmed.fastq in2=short_reads_2-trimmed.fastq out=short_reads_1-trimmed-100x.fastq out2= short_reads_1-trimmed-100x.fastq target=100 mindepth=6 ecc=f 
```

d) Aligning normalized Illumina reads with the Colonization Island sequence

```
# Create a database

bowtie2-build -f ColonizationIsland-LpWF-Masked.fasta dbname

# Align the Paired Short-Reads

bowtie2 -x dbname -1 short_reads_1-trimmed-100x.fastq -2 short_reads_2-trimmed-100x.fastq -S ShortReadsAln-100x.sam --no-unal

#Convert SAM to BAM for sorting

samtools view -S -b ShortReadsAln-100x.sam > ShortReadsAln-100x.bam

#Sort BAM for SNP calling

samtools sort ShortReadsAln-100x.bam ShortReadsAln-sorted-100x.bam
```

e) Evaluating the quality of mapping with Picard

```
java -jar picard.jar CollectAlignmentSummaryMetrics \
	REFERENCE=my_data/ColonizationIsland-LpWF-Masked.fasta \
	INPUT=my_data/ShortReadsAln-sorted-100x.bam \
	OUTPUT=results/Metrics.txt
```


## **7. *In silico* detection of circular and linear contigs.**

a) Split the genome assembly into single contigs 

```
awk '/^>/ {OUT=substr($0,2) ".fasta"}; {print >> OUT; close(OUT)}' HiFiPacBio-circlator.fasta
```

b) Extracting 500 bp from both ends of the contig

```
perl ExtractingEndSeq.pl Contig.fna 500
```

c) Merge both ends of the contig

```
cat contig.rightseq.fa contig.leftseq.fa > Contig-BothSides.fa

# Remove the ">"
```

d) Aligning normalized Illumina reads with the file that has the sequence of both ends of the contig

```
conda activate MappingIlluminaReads

# Create a database

bowtie2-build -f Contig-BothSides.fa dbname-LinearContig

# Align the Paired Short-Reads

bowtie2 -x dbname-LinearContig -1 short_reads_1-trimmed-100x.fastq -2 short_reads_2-trimmed-100x.fastq -S ShortReadsAln-100x-LinearContig.sam --no-unal

#Convert SAM to BAM for sorting

samtools view -S -b ShortReadsAln-100x-LinearContig.sam > ShortReadsAln-100x-LinearContig.bam

#Sort BAM for SNP calling

samtools sort ShortReadsAln-100x-LinearContig.bam ShortReadsAln-sorted-100x-LinearContig.bam
```

## **8. Genome mining for colonization islands in *L. plantarum* genomes.**

a) Using raw data obtained from SRA-NCBI

```
# Installing the SRA Toolkit, Bowtie2, Trim-galore, Minimap2, and htseq-count

conda create --name SRA-DB
conda activate SRA-DB

conda install -c bioconda sra-tools
conda install -c "bioconda/label/cf201901" sra-tools

conda install -c bioconda bowtie2
conda install -c "bioconda/label/broken" bowtie2
conda install -c "bioconda/label/cf201901" bowtie2

conda install -c bioconda trim-galore
conda install -c "bioconda/label/cf201901" trim-galore

conda install -c bioconda htseq
conda install -c "bioconda/label/broken" htseq
conda install -c "bioconda/label/cf201901" htseq

conda install -c bioconda minimap2
conda install -c "bioconda/label/cf201901" minimap2

# Download the raw-reads

prefetch --option-file SRA-DB.txt

# To get the fastq files

fastq-dump -I --split-files SRR390728

# To get the paired-end raw reads

nohup fastq-dump -I --split-files *.sra &

# Trimming the raw data using Trim-galore

./TrimmingSingleEnd-SRA.sh

./TrimmingEndPair-SRA.sh

# # Create a database

bowtie2-build -f ColonizationIsland-LpWF-Masked.fasta dbname

# Align the raw reads vs. the colonization island (ColonizationIsland-LpWF-Masked.fasta)

./AlnSingleEnd-SRA.sh

./AlnPairedEnd-SRA.sh

./AlnSingleEnd-Nano.sh

./AlnSingleEnd-PacBio.sh

# Quantify the reads mapped

./ReadCount.sh
```


## **9. Prediction of the recombination sites.**

a) Installing Yass program, a genomic similarity search tool

```
Please visit: https://github.com/laurentnoe/yass 
```

b) Run Yass to get the Dot plot using the colonization island sequence from LpWF 

```
yass ColonizationIslandLpWF.fna  ColonizationIslandLpWF.mfa    -C 2,-2,-3   -G -5,-2   -E 1e-20   -o LpIsland-output.yop
yass2dotplot.php       LpIsland-output.yop  filename1=""  filename2="" ; open LpIsland-DotPlot.png
```


## **10. Transposable elements annotation and classification into families.**

a) Installing Seqkit and Blast, and python 3

```
conda create --name TEs python=3
conda activate TEs

conda install -c bioconda blast
conda install -c "bioconda/label/cf201901" blast

conda install -c bioconda seqkit
conda install -c "bioconda/label/cf201901" seqkit

```

b) Extracting the TEs amino acid sequences from the fasta file

```
# Download the TEs database from ISfinder (https://github.com/orangeSi/ISfinder_database)

python3 ISfiner.step1.py  step1
python3 ISfiner.step2.py step1.IS.ID.list outdir
rm IS.database.tmp.fa
find outdir/ -type f -name 'IS*.seq.fa'|xargs -L 1 -I {} cat {} >> IS.database.tmp.fa
python3 check.empty.py IS.database.tmp.fa IS.database.fa

# create a list of the TEs IDs saved as TEs-IDs.txt and run:

seqkit grep -n -f TEs-IDs.txt ColonizationIsland-LpWF.fasta > TEs.fasta

# Split the multifasta file in several fasta files

cat TEs.fasta | awk '{if (substr($0, 1, 1)==">") { filename=(substr($0,2) ".fasta")} print $0 > filename}'
```

c) BlastP for TEs annotation based on sequence homology

```
makeblastdb -in ISfinderDB.faa -dbtype prot -out Database.db
blastp -db Database.db -query TEs_aminoseq.fasta -outfmt 6 -evalue 0.001 -num_threads 16 -max_target_seqs 5 -out TEs_aminoseq.blast
```


## **11. Colonization island in other bacteria genera.**

a) Installing HMMER, BLAST, Muscle, Gblocks, biopython, and quicktree 

```
conda create --name aSecTree
conda activate aSecTree

conda install -c bioconda hmmer
conda install -c "bioconda/label/cf201901" hmmer

conda install -c bioconda muscle
conda install -c "bioconda/label/cf201901" muscle

conda install -c bioconda gblocks
conda install -c "bioconda/label/cf201901" gblocks

conda install -c conda-forge biopython
conda install -c "conda-forge/label/cf201901" biopython
conda install -c "conda-forge/label/cf202003" biopython
conda install -c "conda-forge/label/gcc7" biopython

conda install -c bioconda quicktree
conda install -c "bioconda/label/cf201901" quicktree

```

b) Using the BlastP function in PATRIC, to find homologous aSec proteins: SecA2, SecY2, GtfA, GftB, Asp1, Asp2, and Asp3 from LpWF.

```
p3-submit-BLAST \
--max-hits 100000 \
--in-fasta-file aSecLpWF \
--evalue-cutoff 1e-20 \
--workspace-upload-path /userPATRIC@patricbrc.org/home/test \
--db-database BV-BRC \
/userPATRIC@patricbrc.org/home/test test

```

c) Select the best hit per genera: Lower number of contigs, bigger genome size

```
python SelectOneHit-loop.py
```

d) Annotate the selected bacterial genomes using RAST and download them in *.faa

e) Create a genome assemblies database and run HMMER against SecA2 and SecY2 proteins from LpWF 

```
cat *.faa > all_SelectedGenomes.faa
./phmmer --tblout HMMER-SelectedGenomes-SrpA.txt SrpA.fasta all_SelectedGenomes.faa
./phmmer --tblout HMMER-SelectedGenomes-SrpB.txt SrpB.fasta all_SelectedGenomes.faa
```

f) Keep only the column 1 (RAST-IDs)

```
$cut -f 1 HMMER-SelectedGenomes-SrpA.txt > SrpA-Ids
$cut -f 1 HMMER-SelectedGenomes-SrpY.txt > SrpY-Ids
```

g) Extract ~100 proteins upstream and downstream from the HMMER hit

```
./Extract100Seq.sh

The output file "results.txt" will be used to find the aSec proteins using BlastP
```

h) Finding Homologous aSec proteins using BlastP

```
# Make a blast database

perl BlastFormat.pl

# Run BlastP using the aSec proteins as query from LpWF (Asp1.fasta, Asp2.fasta, Asp3.fasta, SecY2.fasta, SecA2.fasta, GftA.fasta, and GftB.fasta)

for i in $(ls *.fasta); do perl BlastP.pl $i $i.txt 0.000001 100;done

# Change the name of the headers for the RAST-Names, using the RAST.Ids file.

perl RASTID_Names.pl aSecProteins-out.faa > aSecProteins-outNames.fasta
```

j) Checking if there are duplicate hits or missing hits for each aSec protein

```
# Create a DocumentA.txt with the RAST.Ids (First column)

cut -f 1 RAST.Ids > DocumentA.txt

# Create a document with a list of missing IDs and duplicate IDs for each Blast-output file (i.e., aSecProteins-outNames.fasta)

python DuplicatesAndMissingIDs.py

# Remove the genomes that do not have a complete aSec protein from all fasta files, creating a list from the Selected Genomes with all the aSec proteins ==> list.txt

seqkit grep -n -f list.txt aSecProteins-outNames.fasta > aSec-Filter.fasta.txt

# Align the resulting fasta files 

for i in $(ls *.fasta.txt); do muscle -in $i -out $i.aln;done

# Trimming the alignments

gblocks /path/to/input.fasta -t=d -e=".gb" -b4=5 -b5=a

- The resulting files are asp1-gb.fasta, asp2-gb.fasta, asp3-gb.fasta, gftA-gb.fasta, gftB-gb.fasta, secA2-gb.fasta, secY2-gb.fasta
```

k) Constructing the final matrix with the aSec proteins to construct the aSec tree.

```
# Extracting the proteins from the asp1-gb.fasta, asp2-gb.fasta, asp3-gb.fasta, gftA-gb.fasta, gftB-gb.fasta, secA2-gb.fasta, secY2-gb.fasta, concatenate them in the same order, and get a final fasta file as output

python ExtractOrderSec.py
```

l) Constructing the aSec tree

```
# Convert the aSec matrix (salida.fasta) to Stockholm format

python stochkolm.py

# Run quicktree

quicktree -in a -out t -boot 1000 aSec.stockholm > aSec.tree

# Visualize in Figtree

(http://tree.bio.ed.ac.uk/software/figtree/)

```

m) Visualizing the genomic islands 

```
The genomic context of the genomic islands was visualized using GeneSpy,
please visit: https://lbbe-dmz.univ-lyon1.fr/GeneSpy/
```

## **11. Co-phylogeny plot.**

a) Installing the programs Quicktree, ete3, Muscle, Gblocks, and Gotree

```
conda create --name Cophylogeny
conda activate Cophylogeny

conda install -c bioconda gotree

conda install -c bioconda quicktree
conda install -c "bioconda/label/cf201901" quicktree

conda install -c conda-forge ete3
conda install -c "conda-forge/label/cf201901" ete3
conda install -c "conda-forge/label/cf202003" ete3

conda install -c bioconda muscle
conda install -c "bioconda/label/cf201901" muscle

conda install -c bioconda gblocks
conda install -c "bioconda/label/cf201901" gblocks

```

b) Aligning the core proteins

```
muscle -in CoreGenome.fasta -out CoreGenome.aln
```

c) Getting the consensus blocks 

```
gblocks CoreGenome.aln -t=d -e=".gb" -b4=5 -b5=a
```

d) Convert the CoreGenome.aln to Stockholm format

```
python stochkolm.py
```

e) Constructing the Core Genome Tree

```
quicktree -in a -out t -boot 1000 CoreGenome.stockholm > CoreGenome.tree
```

f) Midpoint root Core Genome tree

```
python midpoint-root.py CoreGenome.tree > CoreGenome-midpoint.tree
```

g) Collapsing nodes with support <70 in aSec-midpoint.tree and CoreGenome-midpoint.tree

```
gotree collapse support -i aSec-midpoint.tree -s 70 -o aSec-midpoint-collapsed.tree
gotree collapse support -i CoreGenome-midpoint.tree -s 70 -o CoreGenome-midpoint-collapsed.tree
```

h) Construct the co-phylogeny plot and calculate the ParaFit index

```
Run the script "Co-phylogeny plot" in R
```

## **12. SRRPs similarity network.**


a) Run phmmer online to find Orthologs genes using the Bacterial Ensembl Genomes Database

```
Please Visit https://www.ebi.ac.uk/Tools/hmmer/search/phmmer
# Use the Significance E-value -20
# Download the HMMER hits in fasta format
```

b) Construct the SRRPs similarity network

```
EFI-Enzyme similarity tool 
Please visit https://efi.igb.illinois.edu/efi-est/
```

c) Visualize the network using Cytoscape













