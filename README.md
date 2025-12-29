# RNA-Seq_Analysis

## INTRODUCTION
RNA sequencing is a widely used method to study gene expression by measuring RNA levels in cells. It allows us to compare how genes are turned on or off under different biological conditions. RNA-seq is commonly used in cancer research to understand how cells respond to stress, treatment, or changes in their environment.

Prostate cancer usually begins as an androgen-dependent disease, meaning that cancer cells rely on androgen receptor (AR) signaling for growth. Over time, especially under stress or therapy, prostate cancer cells can adapt and become more aggressive and androgen-independent. One such aggressive form is neuroendocrine prostate cancer (NEPC), which is associated with poor prognosis and limited treatment options.

The research study “ONECUT2 is a driver of neuroendocrine prostate cancer” identified the transcription factor ONECUT2 as an important regulator in this process. ONECUT2 suppresses androgen-related gene expression and activates stress-response and hypoxia-related pathways, helping prostate cancer cells survive under unfavorable conditions. Hypoxia, which refers to low oxygen availability, is a common feature of solid tumors and acts as a strong cellular stress that can alter gene expression patterns.

The main aim of this project is to learn and understand the complete RNA-seq analysis workflow, with a particular focus on differential gene expression analysis. For this purpose, RNA-seq data from the LNCaP prostate cancer cell line were used as the primary model system. RNA-seq data from the PC3 cell line were also used during the initial exploratory steps, such as quality control and principal component analysis (PCA), to gain experience in handling multi-sample datasets and to understand overall variation across samples. However, differential gene expression analysis was performed only on LNCaP samples, as using a single cell line is sufficient to learn the RNA-seq workflow and allows a straightforward comparison between normoxia and hypoxia conditions.
___

## DATASET DESCRIPTION

The RNA-seq data used in this project were obtained from the Gene Expression Omnibus (GEO) under the accession GSE106305. This dataset was generated as part of the study “ONECUT2 is a driver of neuroendocrine prostate cancer” and includes RNA-seq experiments performed on human prostate cancer cell lines under different oxygen conditions.

The dataset contains RNA-seq samples from two prostate cancer cell lines:

1. LNCaP (androgen-sensitive)

2. PC3 (androgen-independent)

Cells were cultured under:

1. Normoxia (normal oxygen conditions)

2. Hypoxia (low oxygen conditions)

**LNCaP RNA-seq Samples**
| Condition |	Replicate	|SRA Runs (SRR) |
| ---- | ---- | ---- |
| Normoxia |	Replicate 1	|	SRR7179504, SRR7179505, SRR7179506, SRR7179507 |
| Normoxia |	Replicate 2	|	SRR7179508, SRR7179509, SRR7179510, SRR7179511 |
| Hypoxia |	Replicate 1	| SRR7179520, SRR7179521, SRR7179522, SRR7179523 |
| Hypoxia |	Replicate 2	| SRR7179524, SRR7179525, SRR7179526, SRR7179527 |

**PC3 RNA-seq Samples**
| Condition |	Replicate	|SRA Runs (SRR) |
| ---- | ---- | ---- |
| Normoxia |	Replicate 1	|	SRR7179536 |
| Normoxia |	Replicate 2	|	SRR7179537 |
| Hypoxia |	Replicate 1	|	SRR7179540 |
| Hypoxia |	Replicate 2	|	SRR7179541 |

1. RNA-seq data from both LNCaP and PC3 cell lines were used for initial exploratory analysis, including quality control and PCA.

2. Differential gene expression analysis was performed only on LNCaP samples, comparing normoxia and hypoxia conditions.

3. Multiple SRA runs corresponding to the same biological replicate were merged prior to downstream analysis.

4. All sequencing data are single-end RNA-seq reads.
___

## RNA-SEQ WOEKFLOW

### 1. Downloaded raw sequencing data from SRA

The RNA-seq data for this project were obtained from NCBI SRA. Each sample in GEO is linked to one or more SRR accessions, which store raw sequencing data in .sra format. These .sra files are not directly usable for analysis and must first be downloaded and converted into FASTQ format. 

`prefetch SRR7179504`  This downloads the .sra file for one sequencing run. Repeated this for all SRR accessions listed in the dataset.


```
ls 1_Raw_reads/
  SRR7179504  SRR7179506  SRR7179508  SRR7179510  SRR7179520  SRR7179522  SRR7179524  SRR7179526  SRR7179536  SRR7179540
  SRR7179505  SRR7179507  SRR7179509  SRR7179511  SRR7179521  SRR7179523  SRR7179525  SRR7179527  SRR7179537  SRR7179541`
```


### 2. Converted SRA files to FASTQ format

SRA files are converted into compressed FASTQ files (.fastq.gz) using `fastq-dump.` Each SRR produces a FASTQ file containing sequencing reads.

```

fastq-dump --outdir fastq --gzip --skip-technical --readids \
--read-filter pass --dumpbase --split-3 --clip SRR7179504/SRR7179504.sra`
```
Instead of running these commands repeatedly, the conversion was automated using a Python script `fastq_download.py` that loops over all SRR IDs.

<img width="1919" height="258" alt="Screenshot 2025-12-13 150858" src="https://github.com/user-attachments/assets/a1f3ad65-60d0-4291-8027-5746b9a1441d" />



### 3. Merged technical replicates into sample-wise FASTQ files

Multiple SRR FASTQ files corresponding to one biological replicate are concatenated into a single FASTQ file.

```
cat SRR7179504_pass.fastq.gz SRR7179505_pass.fastq.gz SRR7179506_pass.fastq.gz SRR7179507_pass.fastq.gz  > LNCAP_Normoxia_S1.fastq.gz
cat SRR7179508_pass.fastq.gz SRR7179509_pass.fastq.gz SRR7179510_pass.fastq.gz SRR7179511_pass.fastq.gz  > LNCAP_Normoxia_S2.fastq.gz
cat SRR7179520_pass.fastq.gz SRR7179521_pass.fastq.gz SRR7179522_pass.fastq.gz SRR7179523_pass.fastq.gz  > LNCAP_Hypoxia_S1.fastq.gz
cat SRR7179524_pass.fastq.gz SRR7179525_pass.fastq.gz SRR7179526_pass.fastq.gz SRR7179527_pass.fastq.gz  > LNCAP_Hypoxia_S2.fastq.gz

mv SRR7179536_pass.fastq.gz PC3_Normoxia_S1.fastq.gz
mv SRR7179537_pass.fastq.gz PC3_Normoxia_S2.fastq.gz
mv SRR7179540_pass.fastq.gz PC3_Hypoxia_S1.fastq.gz
mv SRR7179541_pass.fastq.gz PC3_Hypoxia_S2.fastq.gz

```

### 4. Initial quality check (FastQC)

FastQC utilised on the sample-wise FASTQ files to assess sequencing quality. FastQC checks various parameters like base quality scores, GC content, adapter contamination, overrepresented sequences. 

Command
```
fastqc 3_Fastq_samplewise/*.fastq.gz \
-o 4_Fastqc_results/ --threads 8
```
<img width="1919" height="469" alt="Screenshot 2025-12-13 190450" src="https://github.com/user-attachments/assets/09a8a039-25a0-4bc6-9b52-91178ed0cee3" />

I have attached one file, named `LNCAP_Hypoxia_S1.fastq.gz FastQC Report.pdf`, you can see there the QC result of that sample.

MultiQC is used to collect and summarize all FastQC results into a single consolidated report.
```
multiqc 4_Fastqc_results/ -o 4_Fastqc_results/
```

### 5. Read Trimming

Raw RNA-seq reads often contain low-quality bases at the ends of reads, Adapter contamination, Sequencing artifacts. Read trimming removes these unwanted portions of reads to improve data quality before alignment. Trimming was performed using Trimmomatic on single-end reads.

```
java -jar Software/Trimmomatic-0.39/trimmomatic-0.39.jar SE -threads 4 -phred33 \
  fastq/SRR7179504_pass.fastq.gz \
  fastq/SRR7179504_trimmed.fastq.gz \
  TRAILING:10
```
For each input FASTQ file, Trimmomatic produces a trimmed FASTQ file. Only trimmed reads are used for all downstream steps. 

### 6. Genome Alignment Using HISAT2

After trimming, the RNA-seq reads are aligned to the human reference genome (GRCh38). Alignment means determining where each sequencing read originates in the genome. This step converts raw sequencing reads into genome-mapped reads, which are required for gene-level quantification. 

Downloaded HISAT2 prebuilt GRCh38 genome index
```
wget https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz
tar -xvzf grch38_genome.tar.gz
```

Downloaded Ensembl GTF annotation:
```
wget https://ftp.ensembl.org/pub/release-114/gtf/homo_sapiens/Homo_sapiens.GRCh38.114.gtf.gz
gunzip Homo_sapiens.GRCh38.114.gtf.gz

```
A GTF (Gene Transfer Format) file contains information about Gene locations, Exon boundaries, Transcript structures, Gene identifiers. RNA-seq reads align to the genome, not directly to genes. To count how many reads belong to each gene, we need gene coordinates. The GTF file provides these coordinates. Without a GTF file gene-level quantification is not possible.

For one sample:
```
hisat2 -q -x grch38/genome -U fastq/LNCAP_Hypoxia_S1_trimmed.fastq.gz | \
  samtools sort -o alignedreads/LNCAP_Hypoxia_S1.bam

samtools index alignedreads/LNCAP_Hypoxia_S1.bam
```
To avoid running alignment commands manually for each sample, alignment was automated using a shell script. Script used: `hisat2alignment.sh`. 
It will loops over all trimmed FASTQ files, aligns each sample using HISAT2, sorts BAM files using Samtools, indexes BAM files, logs runtime information.

<img width="1888" height="140" alt="Screenshot 2025-12-29 200130" src="https://github.com/user-attachments/assets/689862da-717c-47ff-b172-a2afae43c6c5" />


### 7. Assessing quality of aligned reads

After aligning RNA-seq reads to the reference genome, it is important to check how good the alignment actually is. This step evaluates whether reads mapped correctly, evenly, and in a biologically reasonable way. Qualimap is used to assess the quality of RNA-seq alignments stored in BAM files.

```
qualimap rnaseq -bam alignedreads/LNCAP_Hypoxia_S1.bam \
  -gtf Homo_sapiens.GRCh38.114.gtf \
  -outdir rnaseq_qc_results --java-mem-size=8G
```

For each sample, Qualimap generates a qualimapReport.html, an interactive QC report. I have attached one file, named `LNCAP_Hypoxia_S1.fastq.gz FastQC Report.pdf`, you can see there the QC result of that sample.

