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

The next step in the RNA-seq workflow is the differential expression analysis. The goal of differential expression testing is to determine which genes are expressed at different levels between conditions. These genes can offer biological insight into the processes affected by the conditions of interest. DESeq2 is a widely used Bioconductor package in R designed for differential gene expression analysis of high-throughput sequencing data, such as RNA-seq. So, required packages installed.

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("DESeq2", "apeglm", "EnhancedVolcano", "pheatmap", "RColorBrewer"))
install.packages(c("tidyverse", "ggrepel"))
```

#### 1. Loading required libraries

```r
library(DESeq2)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(matrixStats)
library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)
library(fgsea)
```
DESeq2 - Core package used for normalization and differential gene expression analysis of RNA-seq count data.

dplyr, tidyverse,  - Used for efficient data manipulation, filtering, formatting, and handling sample and gene information.

ggplot2, pheatmap, RColorBrewer - Used for visualization of results, including PCA plots, volcano plots, and heatmaps.

clusterProfiler, org.Hs.eg.db, ReactomePA, fgsea - Used for functional annotation and pathway enrichment analysis of differentially expressed genes.

 #### 2. Loading Count Matrix 
```r
#Read count matrix

counts <- read.csv("GSE106305_counts_matrix.csv",
row.names = "Geneid",
stringsAsFactors = FALSE)


counts <- counts[, sort(colnames(counts))]
#View the first few rows
head(counts)

#Rows = genes, Columns = RNA-seq samples, Values = raw read counts (NOT normalized)
```


#### 3. Creating Sample Metadata

```{r}
condition <- factor(c(
rep("LNCAP_Hypoxia", 2),
rep("LNCAP_Normoxia", 2),
rep("PC3_Hypoxia", 2),
rep("PC3_Normoxia", 2)
))


colData <- data.frame(condition)
rownames(colData) <- colnames(counts)
head(colData)

#Metadata encodes the experimental design. DESeq2 estimates expression changes relative to these biological conditions.
```


#### 4. Creating DESeq2 Dataset and Inspect Zero-Count Genes

```r
dds <- DESeqDataSetFromMatrix(
countData = counts,
colData = colData,
design = ~ condition
)
#This object links gene counts with experimental design, enabling statistical inference of differential expression.


keep <- rowSums(counts(dds) >= 10) >= 2
#Keep genes with at least 10 reads in at least 2 samples

dds <- dds[keep, ]
```
Genes expressed in almost no samples carry little biological signal. Removing them improves dispersion estimation. Avoids false positives driven by noise.


#### 5. Differential expression analysis

```{r}
dds <- DESeq(dds)
#Here DESeq2 Estimates size factors (library depth), Estimates gene-wise dispersion, Fits negative binomial GLMs, Performs hypothesis testing


#Normalized counts (used for visualization only)
norm_counts <- counts(dds, normalized = TRUE)
write.csv(norm_counts, "Normalized_counts.csv")

#Never use normalized counts for DE testing. DESeq2 already handled normalization internally.
```



#### 6. Exploratory analysis: PCA
```{r}

#Variance-stabilizing transformation
vsd <- vst(dds, blind = TRUE)
#Raw counts are highly skewed. VST makes variance roughly constant across expression levels with preserving biological differences.
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))


p_pca <- ggplot(pca_data, aes(PC1, PC2, color = condition)) +
geom_point(size = 3) +
labs(x = paste0("PC1: ", percentVar[1], "%"),
y = paste0("PC2: ", percentVar[2], "%"),
title = "PCA (VST transformed)") +
theme_minimal()


ggsave("PCA_plot.png", p_pca, width = 6, height = 5, dpi = 300)

print(p_pca)
dev.off()
```
PCA answers "what dominates variation?"
# pca plot

#### 7. Sample-to-sample distance heatmap
```{r}
sample_dists <- dist(t(assay(vsd)))
dist_matrix <- as.matrix(sample_dists)


png("sample_distance_heatmap.png", width = 1200, height = 1000, res = 300)
p <- pheatmap(dist_matrix,
clustering_distance_rows = sample_dists,
clustering_distance_cols = sample_dists,
col = colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
main = "Sample-to-sample distances")

print(p)
dev.off()
```
Dark blue = similar transcriptomes

Light color = dissimilar transcriptomes
# plot

#### 8. Expression distribution diagnostics
```{r}
png("density_raw_vs_vst.png", width = 2000, height = 2000, res = 300)
par(mfrow = c(2, 2))


plot(density(counts(dds)[,1]), main = "Raw counts", xlab = "Expression")
plot(density(assay(vsd)[,1]), main = "VST counts", xlab = "Expression")


dev.off()

```
Raw counts:

heavy right tail

dominated by few highly expressed genes

VST counts:

approximately normal

comparable across samples

This confirms transformation behaved correctly.
# plot


#### 9. Highly variable gene heatmap
```{r}
rv <- rowVars(assay(vsd))
top_genes <- order(rv, decreasing = TRUE)[1:40]


mat <- assay(vsd)[top_genes, ]
mat <- mat - rowMeans(mat)


png("top_variable_genes_heatmap.png", width = 1200, height = 1200, res = 300)
pheatmap(mat,
color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(255),
fontsize_row = 6,
main = "Top 40 most variable genes")
dev.off()
```
These genes:

drive major differences between samples

often include stress-response, cell-cycle, or condition-specific genes

This heatmap helps confirm:

consistency of biological signal

absence of sample-level artifacts
# plot

```{r}
dds_lncap <- dds[, grepl("LNCAP", colnames(dds))]
dds_lncap$condition <- droplevels(dds_lncap$condition)
dds_lncap$condition <- relevel(dds_lncap$condition, ref = "LNCAP_Normoxia")

```


#### 10. Extracted differential expression results
```{r}
dds_lncap <- DESeq(dds_lncap)


res_lncap <- results(dds_lncap,
contrast = c("condition",
"LNCAP_Hypoxia",
"LNCAP_Normoxia"))
write.csv(as.data.frame(res_lncap), "DEGs_LNCAP.csv")

```

#### 11. MA plot

```{r}

png("MA_plot_LNCAP.png", width = 800, height = 600, res = 150)
plotMA(res_lncap, ylim = c(-5, 5))
dev.off()
```
# plot

#### 12. Volcano plot
```{r}

res_df <- as.data.frame(res_lncap) %>%
na.omit() %>%
mutate(regulation = case_when(
padj < 0.05 & log2FoldChange > 1 ~ "Upregulated",
padj < 0.05 & log2FoldChange < -1 ~ "Downregulated",
TRUE ~ "Not significant"
))


p_volcano <- ggplot(res_df, aes(log2FoldChange, -log10(padj), color = regulation)) +
geom_point(alpha = 0.6) +
theme_minimal() +
labs(title = "Volcano plot: LNCaP hypoxia")


ggsave("Volcano_plot_LNCAP.png", p_volcano, width = 6, height = 5, dpi = 300)
```
# plot

### 10. Pathway Analysis

RNA-seq differential expression produces lists of genes. However, biology does not operate gene-by-gene.
Instead, it operates through pathways, signaling cascades, coordinated gene programs.
Pathway analysis helps to answer, Which biological processes or programs are altered between conditions?

There are different methods utilised to answer this question. Here we are analysing, using three different approaches. 
All three methods aim to answer the same high-level question:
Which biological pathways are affected between conditions?
But they do so using different assumptions and levels of information.

| Method	| Core idea |	Information used |
| ---------	| -------- |	---------- |
| GSEA	| Distribution of genes |	All genes (ranked) |
| ORA	| Counting genes	| Only significant genes |
| fgsea + Hallmark |	Fast, robust GSEA	| All genes, curated sets |


#### PART 1: Pre-ranked GSEA (Reactome)

##### Step 1: Load required libraries

```r
library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)
library(dplyr)
```



##### Step 2: Convert gene IDs (ENSEMBL → ENTREZ)

**Why:** Reactome pathways are indexed by ENTREZ IDs

```r
res_lncap$ENSEMBL <- rownames(res_lncap)

id_map <- bitr(
  res_lncap$ENSEMBL,
  fromType = "ENSEMBL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)
```


##### Step 3: Merge mapping and remove duplicates

Each gene must contribute once; duplicates distort enrichment.

```r
res_mapped <- res_lncap %>%
  left_join(id_map, by = "ENSEMBL") %>%
  filter(!is.na(ENTREZID)) %>%
  distinct(ENTREZID, .keep_all = TRUE)
```

---

##### Step 4: Create ranked gene list 

Ranking metric: `log2FoldChange`

```r
gene_ranking <- res_mapped$log2FoldChange
names(gene_ranking) <- res_mapped$ENTREZID

gene_ranking <- sort(gene_ranking, decreasing = TRUE)
```

*  Top genes are genes induced by condition
* Bottom genes are genes repressed by condition

---

##### Step 5: Run Reactome GSEA

```r
gsea_reactome <- gsePathway(
  geneList     = gene_ranking,
  organism     = "human",
  pvalueCutoff = 0.05,
  verbose      = FALSE
)
```

---

## Step 6: Inspect enriched pathways

```r
head(gsea_reactome@result)
```

**Key interpretation rules:**

* NES > 0 → pathway activated
* NES < 0 → pathway suppressed

---

#### PART 2: Over-Representation Analysis (ORA)

* This particular analysis talks about, Are significant DEGs over-represented in known pathways?


##### Step 1: Define significant genes


```r
sig_genes <- res_mapped %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1) %>%
  pull(ENTREZID)
```



##### Step 2: Run ORA using Reactome

```r
ora_reactome <- enrichPathway(
  gene         = sig_genes,
  organism     = "human",
  pvalueCutoff = 0.05
)
```


##### Step 3: Visualize ORA results

```r
dotplot(ora_reactome, showCategory = 20)
```

ORA ignores genes just below the cutoff, depends heavily on arbitrary thresholds, misses subtle
coordinated shifts. Therefore, ORA should be treated as a supporting or confirmatory analysis, not the primary one.
# plot
---
#### PART 3: fgsea + Hallmark 

fgsea is a fast and efficient implementation of the GSEA algorithm. It asks the same question as GSEA but
scales better and is widely used. 
Hallmark pathways are summarize core biological programs and manually curated. 
It's reduce redundancy and highlight dominant biology.
This makes fgsea + Hallmark ideal for final biological conclusions.


##### Step 1: Load fgsea and Hallmark gene sets

Downloaded GMT file manually (h.all.v7.0.symbols.gmt.txt), which contains the Hallmark gene set collection from MSigDB version 7.0. https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp

```r
library(fgsea)

hallmark_sets <- gmtPathways("h.all.v7.0.symbols.gmt.txt")n```
```



##### Step 2: Prepared ranked list

```r
# Convert ENSEMBL to SYMBOL
symbol_map <- bitr(
  res_lncap$ENSEMBL,
  fromType = "ENSEMBL",
  toType   = "SYMBOL",
  OrgDb    = org.Hs.eg.db
)

res_symbol <- res_lncap %>%
  left_join(symbol_map, by = c("ENSEMBL" = "ENSEMBL")) %>%
  filter(!is.na(SYMBOL)) %>%
  distinct(SYMBOL, .keep_all = TRUE)

ranked_symbols <- res_symbol$log2FoldChange
names(ranked_symbols) <- res_symbol$SYMBOL

ranked_symbols <- sort(ranked_symbols, decreasing = TRUE)
````


##### Step 3: Run fgsea

```r
fgsea_res <- fgsea(
  pathways = hallmark_sets,
  stats    = ranked_symbols,
  minSize  = 15,
  maxSize  = 500,
  nperm    = 1000
)
```



##### Step 4: Inspect top Hallmark pathways

```r
fgsea_res %>%
  arrange(padj) %>%
  select(pathway, NES, padj) %>%
  head()
```



##### Step 5: Visualized Hallmark enrichment

```r
library(ggplot2)
library(stringr)

fgsea_res %>%
  arrange(padj) %>%
  slice(1:15) %>%
  mutate(pathway = str_remove(pathway, "HALLMARK_")) %>%
  ggplot(aes(x = NES, y = reorder(pathway, NES), fill = padj < 0.05)) +
  geom_col() +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "Hallmark pathways altered",
    x = "Normalized Enrichment Score",
    y = NULL
  )
```

# plot
---


Using a combination of pre-ranked GSEA, ORA, and fgsea with Hallmark gene sets, we analysed the transcriptional consequences of hypoxia in LNCaP prostate cancer cells. Pathway-level analysis comparing LNCaP cells grown under normoxia and hypoxia reveals a clear  transcriptional response to low-oxygen stress. Using fgsea with Hallmark gene sets, we observed strong positive enrichment of pathways related to hypoxia response, glycolysis, and cellular stress adaptation, indicating coordinated upregulation of genes involved in oxygen sensing and metabolic reprogramming. In contrast, pathways associated with differentiated cellular functions show relative suppression under hypoxic conditions.

Overall, fgsea with Hallmark gene sets provides a robust and interpretable summary of hypoxia-induced biological programs in LNCaP cells.



