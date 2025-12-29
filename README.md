# RNA-Seq_Analysis

<ins>**INTRODUCTION**</ins>

RNA sequencing is a widely used method to study gene expression by measuring RNA levels in cells. It allows us to compare how genes are turned on or off under different biological conditions. RNA-seq is commonly used in cancer research to understand how cells respond to stress, treatment, or changes in their environment.

Prostate cancer usually begins as an androgen-dependent disease, meaning that cancer cells rely on androgen receptor (AR) signaling for growth. Over time, especially under stress or therapy, prostate cancer cells can adapt and become more aggressive and androgen-independent. One such aggressive form is neuroendocrine prostate cancer (NEPC), which is associated with poor prognosis and limited treatment options.

The research study “ONECUT2 is a driver of neuroendocrine prostate cancer” identified the transcription factor ONECUT2 as an important regulator in this process. ONECUT2 suppresses androgen-related gene expression and activates stress-response and hypoxia-related pathways, helping prostate cancer cells survive under unfavorable conditions. Hypoxia, which refers to low oxygen availability, is a common feature of solid tumors and acts as a strong cellular stress that can alter gene expression patterns.

The main aim of this project is to learn and understand the complete RNA-seq analysis workflow, with a particular focus on differential gene expression analysis. For this purpose, RNA-seq data from the LNCaP prostate cancer cell line were used as the primary model system. RNA-seq data from the PC3 cell line were also used during the initial exploratory steps, such as quality control and principal component analysis (PCA), to gain experience in handling multi-sample datasets and to understand overall variation across samples. However, differential gene expression analysis was performed only on LNCaP samples, as using a single cell line is sufficient to learn the RNA-seq workflow and allows a straightforward comparison between normoxia and hypoxia conditions.
___

<ins>**DATASET DESCRIPTION**</ins>

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

<ins>**RNA-SEQ WOEKFLOW**</ins>

1. Downloaded raw sequencing data from SRA

The RNA-seq data for this project were obtained from NCBI SRA. Each sample in GEO is linked to one or more SRR accessions, which store raw sequencing data in .sra format. These .sra files are not directly usable for analysis and must first be downloaded and converted into FASTQ format. 

`prefetch SRR7179504`  This downloads the .sra file for one sequencing run. Repeated this for all SRR accessions listed in the dataset.


```
ls 1_Raw_reads/
  SRR7179504  SRR7179506  SRR7179508  SRR7179510  SRR7179520  SRR7179522  SRR7179524  SRR7179526  SRR7179536  SRR7179540
  SRR7179505  SRR7179507  SRR7179509  SRR7179511  SRR7179521  SRR7179523  SRR7179525  SRR7179527  SRR7179537  SRR7179541`
```


2. Converted SRA files to FASTQ format

SRA files are converted into compressed FASTQ files (.fastq.gz) using `fastq-dump.` Each SRR produces a FASTQ file containing sequencing reads.

```

fastq-dump --outdir fastq --gzip --skip-technical --readids \
--read-filter pass --dumpbase --split-3 --clip SRR7179504/SRR7179504.sra`
```
Instead of running these commands repeatedly, the conversion was automated using a Python script `fastq_download.py` that loops over all SRR IDs.

<img width="1919" height="258" alt="Screenshot 2025-12-13 150858" src="https://github.com/user-attachments/assets/a1f3ad65-60d0-4291-8027-5746b9a1441d" />



3. Merged technical replicates into sample-wise FASTQ files

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

4. Initial quality check (FastQC)

FastQC is run on the sample-wise FASTQ files to assess sequencing quality. FastQC checks various parameters like base quality scores, GC content, adapter contamination, overrepresented sequences. 

Command
```
fastqc 3_Fastq_samplewise/*.fastq.gz \
-o 4_Fastqc_results/ --threads 8
```
