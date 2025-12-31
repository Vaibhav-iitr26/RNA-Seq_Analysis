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

```
prefetch SRR7179504
```
This downloads the .sra file for one sequencing run. Repeated this for all SRR accessions listed in the dataset.


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

For each sample, Qualimap generates a qualimapReport.html, an interactive QC report. I have attached one file, named `Qualimap report_ LNCAP_Hypoxia_S1_RNA Seq QC.pdf`, you can see there the qualimap report.

To apply the same quality checks to all samples, Qualimap was run using a shell script. Script used: `run_qualimap.sh`.


### 8. Gene-level Read Quantification Using featureCounts

After alignment and alignment quality checking, each RNA-seq read has a known position in the genome.
In this step, aligned reads are assigned to genes and counted. The output of this step is a gene × sample count matrix, which is the direct input for differential expression analysis. featureCounts is a fast and widely used read-counting tool. 

For a single sample:
```
featureCounts \
  -a Homo_sapiens.GRCh38.114.gtf \
  -o LNCAP_Hypoxia_S1_featurecounts.txt \
  LNCAP_Hypoxia_S1.bam

```
For each sample, featureCounts generates a gene count file. Because featureCounts must be run on multiple BAM files, read counting was automated using a shell script named `featurecounts.sh`.

After running featureCounts, you obtain one count file per sample. Each file contains gene IDs and the number of reads mapped to that gene for that specific sample. So, all individual count files are need to be merged into a single count matrix, where, Rows are genes, Columns are samples, Values are actual raw read counts. To obtain this matrix, python script `countsmatrix_wholedata.py` used. The final output is a combined count matrix stored in file `GSE106305_counts_matrix.csv`. This combined matrix is the direct input for differential expression analysis.


### 9. Differential Gene Expression

The next step in the RNA-seq workflow is the differential expression analysis. The goal of differential expression testing is to determine which genes are expressed at different levels between conditions. These genes can offer biological insight into the processes affected by the condition(s) of interest. DESeq2 is a widely used Bioconductor package in R designed for differential gene expression analysis of high-throughput sequencing data, such as RNA-seq. So, required packages installed.

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("DESeq2", "apeglm", "EnhancedVolcano", "pheatmap", "RColorBrewer"))
install.packages(c("tidyverse", "ggrepel"))
```

#### 1. Loading required libraries

```r
library(DESeq2)
library(dplyr)
library(tibble)
library(tidyverse)
library(data.table)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(matrixStats)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(fgsea)
library(forcats)
library(stringr)
```
DESeq2 - Core package used for normalization and differential gene expression analysis of RNA-seq count data.

dplyr, tibble, tidyverse, data.table, forcats, stringr - Used for efficient data manipulation, filtering, formatting, and handling sample and gene information.

ggplot2, ggrepel, pheatmap, RColorBrewer, matrixStats - Used for visualization of results, including PCA plots, volcano plots, and heatmaps.

clusterProfiler, org.Hs.eg.db, ReactomePA, fgsea - Used for functional annotation and pathway enrichment analysis of differentially expressed genes.

 #### 2. Loading Count Matrix and Creating Sample Metadata
```r
#Loads the combined gene × sample raw count matrix. Gene IDs are used as row names.
raw_counts <- read.csv(
  "8_Read_counts/GSE106305_counts_matrix.csv",
  row.names = "Geneid",
  stringsAsFactors = FALSE
)


#Ensure consistent sample order
raw_counts <- raw_counts[, sort(colnames(raw_counts))]


#Checks data dimensions and total read counts per sample.
dim(raw_counts)
colSums(raw_counts)


#Defines experimental conditions for each sample.
condition <- c(
rep("LNCAP_Hypoxia", 2),
rep("LNCAP_Normoxia", 2),
rep("PC3_Hypoxia", 2),
rep("PC3_Normoxia", 2)
)


#Creates sample metadata aligned with the count matrix.
colData <- data.frame(condition)
rownames(colData) <- colnames(raw_counts)
colData

```

#### 3. Creating DESeq2 Dataset and Inspect Zero-Count Genes

```r
#Creates a DESeq2 dataset by combining the raw count matrix with sample metadata and specifying the experimental design based on the condition variable.
dds <- DESeqDataSetFromMatrix(
countData = raw_counts,
colData   = colData,
design    = ~ condition
)

#Extracts the raw gene-level read counts from the DESeq2 dataset object. 
count_matrix <- counts(dds)

#Calculates the number of samples in which each gene has zero read counts.
zero_counts_per_gene <- rowSums(count_matrix == 0)
table(zero_counts_per_gene)
```


#### 4. 

```r
annotation <- fread("Data/GRCh38annotation.csv")

annotation$Geneid <- sub("\\..*$", "", annotation$Geneid)

counts_df <- raw_counts %>%
rownames_to_column("Geneid") %>%
mutate(Geneid = sub("\\..*$", "", Geneid))

annotated_counts <- left_join(counts_df, annotation, by = "Geneid") %>%
    dplyr::select(Geneid, Genesymbol, Genebiotype, 
           LNCAP_Hypoxia_S1, LNCAP_Hypoxia_S2, LNCAP_Normoxia_S1, LNCAP_Normoxia_S2, 
           PC3_Hypoxia_S1, PC3_Hypoxia_S2, PC3_Normoxia_S1, PC3_Normoxia_S2)

biotypes_to_keep <- c(
"protein_coding",
"IG_J_gene","IG_V_gene","IG_C_gene","IG_D_gene",
"TR_D_gene","TR_C_gene","TR_V_gene","TR_J_gene"
)

filtered_counts <- annotated_counts %>%
filter(Genebiotype %in% biotypes_to_keep)

zero_counts <- rowSums(filtered_counts[, 4:11] == 0)
filtered_counts <- filtered_counts[zero_counts < 7, ]

fwrite(filtered_counts, "filtered_biotype_nozero_count_matrix.csv")

```
[RNA-seq Differential Expression.pdf](https://github.com/user-attachments/files/24397034/RNA-seq.Differential.Expression.pdf)



