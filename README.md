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



## **3. Comparative genomics of *L. plantarum* genome assemblies and their functional annotation based on RAST subsystems.**

```
The comparison was performed using the Function-based comparison tool from RAST
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

conda install -c bioconda bowtie2
conda install -c "bioconda/label/broken" bowtie2
conda install -c "bioconda/label/cf201901" bowtie2

conda install -c bioconda bbmap
conda install -c "bioconda/label/cf201901" bbmap

conda activate MappingIlluminaReads

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


## **8. *L. plantarum* phylogenetic tree.**

a) Installing MrBayes

```
conda create --name LpTree

conda install -c bioconda mrbayes
conda install -c "bioconda/label/cf201901" mrbayes

conda activate LpTree

```

b) High throughput RAST annotation and calculating the core proteome

```
# Please visit: 
Gutiérrez-García, Karina, et al. "Cycad coralloid roots contain bacterial communities including cyanobacteria and Caulobacter spp. that encode niche-specific biosynthetic gene clusters." Genome Biology and Evolution 11.1 (2019): 319-334.

# Use the Lp.Ids file to obtain the core proteome
```

c) Constructing the L. plantarum phylogenetic tree

```
#Convert the LpCoreProteome.fasta in nexus format

Please visit https://doua.prabi.fr/software/seaview

# Add the MrBayes.fa script at the end of the file LpCoreProteome-nexus.nxs

# Run MrBayes
nohup My/path/mrbayes/MrBayes/src/mb -i LpCoreProteome-nexus.nxs > LpCoreProteome-nexus.log&

# Visualize the Lp tree using Iroki

Please visit: https://www.iroki.net/viewer
```


## **8. Genome mining for colonization islands in *L. plantarum* genomes.**

a) Using Lp Genome assemblies 

```
#Installing HMMER

conda create --name LpIsland
conda install -c bioconda hmmer
conda install -c "bioconda/label/cf201901" hmmer
conda activate LpIsland

# Download from RAST the amino acid files for each L. plantarum genome from the Lp.Ids

# Create an Lp database

cat *.faa > all_Lpgenomes_prot.faa
formatdb -i all_genomes_prot.faa -p T -V

# Get all the queries from the colonization island from LpWF

awk '/^>/ {OUT=substr($0,2) ".fasta"}; {print >> OUT; close(OUT)}' ColonizationIslandLpWF.fasta

# Run HMMER

./runHMMER.sh

# Keep only the first and third column

cut -f 1,3 HMMER.out > HMMER-1-3.txt

# Ordering the HMMER-1-3.txt table from ascending to descending

python OrderHMMERTable.py

# Select the lowest e-value per hit

perl HMMER-Evalue.pl

# Adding missing Lp genomes (Genomes that did not have any HMMER- hits) to the file out-HMMER-1-3.order.txt based on the Lp-DB

bash CompletingLpDB.sh

The input files are:
	- Order_LpTree.txt
	- out-HMMER-1-3.order.txt

# Fill empty spaces in the second column with 100

python Adding100.py

# Ordering the table final.results100.txt according to the Lp Tree

./LpTree-Final.sh

The input files are:
	- Order_LpTree.txt
	- final.results100.txt
```

b) Using raw data obtained from SRA-NCBI

```
# Installing the SRA Toolkit, Bowtie2, Trim-galore, Minimap2, and htseq-count

conda create --name SRA-DB
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

conda activate SRA-DB

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

conda install -c bioconda blast
conda install -c "bioconda/label/cf201901" blast

conda install -c bioconda seqkit
conda install -c "bioconda/label/cf201901" seqkit

conda activate TEs
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
conda install -c bioconda hmmer
conda install -c "bioconda/label/cf201901" hmmer

conda install -c bioconda blast
conda install -c "bioconda/label/cf201901" blast

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

conda activate aSecTree
```

b) phmmer to find Orthologs genes using the Bacterial Ensembl Genomes Database

```
-E 0.0000000001 --domE 0.003 --incE 0.0000000001 --incdomE 0.003 --mx BLOSUM62 --pextend 0.4 --popen 0.02 --seqdb ensemblgenomes
```

c) Download the bacterial genomes with HMMER-hits (threshold e-value equal to or less than e-20 and one genome per bacterial genera)

```
Please visit: https://bacteria.ensembl.org/info/data/ftp/index.html
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

# Aling the resulting fasta files 

for i in $(ls *.fasta.txt); do muscle -in $i -out $i.aln;done

# Trimming the alignments

gblocks /path/to/input.fasta -t=d -e=".gb" -b4=5 -b5=a

- The resulting files are asp1-gb.fasta, asp2-gb.fasta, asp3-gb.fasta, gftA-gb.fasta, gftB-gb.fasta, secA2-gb.fasta, secY2-gb.fasta
```

k) Constructing the final matrix with the aSec proteins to construct the aSec tree.

```-
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

## **11. SRRPs similarity network.**

a) Installing HMMER

```
conda create --name SRRPs

conda install -c bioconda hmmer
conda install -c "bioconda/label/cf201901" hmmer

conda activate SRRPs
```

b) phmmer to find Orthologs genes using the Bacterial Ensembl Genomes Database

```
-E 0.0000000001 --domE 0.003 --incE 0.0000000001 --incdomE 0.003 --mx BLOSUM62 --pextend 0.4 --popen 0.02 --seqdb ensemblgenomes
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

SRRPs similarity network. To evaluate the sequence homology among SRRPs from diverse bacterial genera, a similarity network was constructed. Homologs of SrpA and SrpB were identified using pHMMER v3.3.2 (REF) with the Bacterial Ensembl Genomes Database. All the HMMER hits with a threshold e-value equal to or less than e-20 were selected and its adhesin protein sequence was extracted and used to create a SRRPs database. The SRRPs similarity network was constructed using the EFI-Enzyme similarity tool (REF) with the default parameters. The network was visualized using Citoscape v3.10.1 (REF).
















