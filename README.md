# A Workflow for the Detection and Translatome-Guided Prioritization of Cancer Fusion Transcripts

This repository accompanies the book chapter:

> **A Workflow for the Detection and Translatome-Guided Prioritization of Cancer Fusion Transcripts**

The workflow combines RNA-seq-based fusion transcript detection with translatome profiling using matched Ribo-seq and RNA-seq datasets. The aim is to detect candidate fusion transcripts, extract tumor-specific fusion partner genes, and prioritize those partners according to translational regulation patterns.

---

## Workflow overview

```text
RNA-seq FASTQ download
        ↓
Read quality control
        ↓
STAR genome index construction
        ↓
STAR alignment for fusion detection
        ↓
Arriba / STAR-Fusion via nf-core/rnafusion
        ↓
Fusion partner gene extraction
        ↓
Tumor-specific partner filtering
        ↓
Ribo-seq / RNA-seq count processing
        ↓
anota2seq differential translation analysis
        ↓
Translatome-guided fusion partner prioritization
```

---

# 3.1 Fusion transcript detection

## Step 1. Download RNA-seq reads from SRA

RNA-seq reads are downloaded from SRA using **SRA Toolkit v3.4.1** and converted into FASTQ format.

```bash
fasterq-dump --skip-technical -p -e 8 -O ~/GSE112705/fastq SRRXXXXXXX
```

For multiple accessions, the same command can be applied iteratively using a text file containing one SRR accession per line.

```bash
while read SRR
do
  fasterq-dump \
    --skip-technical \
    -p \
    -e 8 \
    -O ~/GSE112705/fastq \
    "$SRR"
done < SRR_Acc_List.txt
```

The file `SRR_Acc_List.txt` should contain one SRR accession per line.

The resulting FASTQ files are then compressed to reduce storage requirements.

```bash
for fq in ~/GSE112705/fastq/*.fastq
do
  gzip "$fq"
done
```

---

## Step 2. Assess raw read quality with FastQC

Raw FASTQ files are assessed using **FastQC v0.12.1**. To evaluate sequencing quality before alignment, run:

```bash
fastqc -t 8 -o ~/GSE112705/fastqc ~/GSE112705/fastq/*.fastq.gz
```

FastQC reports should be inspected for per-base sequence quality, GC-content distribution, sequence duplication, adapter contamination, and overrepresented sequences. Samples with acceptable quality profiles are retained for downstream analysis. If adapter contamination or low-quality terminal bases are observed, trimming should be performed before alignment.

---

## Step 3. Build the STAR genome index

The GRCh38 reference genome and GENCODE v46 annotation are downloaded and used to build a STAR genome index.

```bash
mkdir -p ~/references/GRCh38/gencode_v46/gencode
cd ~/references/GRCh38/gencode_v46/gencode

# Genome FASTA
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh38.primary_assembly.genome.fa.gz

# GTF annotation
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.annotation.gtf.gz

gunzip *.gz

mkdir -p ~/references/GRCh38/gencode_v46/star

STAR \
  --runThreadN 24 \
  --runMode genomeGenerate \
  --genomeDir ~/references/GRCh38/gencode_v46/star \
  --genomeFastaFiles ~/references/GRCh38/gencode_v46/gencode/GRCh38.primary_assembly.genome.fa \
  --sjdbGTFfile ~/references/GRCh38/gencode_v46/gencode/gencode.v46.annotation.gtf \
  --sjdbOverhang 74
```

The `--sjdbOverhang` parameter should be set to the read length minus one. The STAR index is generated once and reused for all samples processed with the same reference and annotation.

---

## Step 4A. Run STAR alignment for Arriba-compatible chimeric detection

```bash
STAR \
  --genomeDir ~/references/GRCh38/gencode_v46/star \
  --readFilesIn ~/GSE112705/fastq/SRRXXXXXXX.fastq.gz \
  --readFilesCommand zcat \
  --runThreadN 8 \
  --outFileNamePrefix ~/GSE112705/star_prep/SRRXXXXXXX/SRRXXXXXXX. \
  --outSAMtype BAM Unsorted \
  --outSAMunmapped Within \
  --outBAMcompression 0 \
  --outFilterMultimapNmax 50 \
  --peOverlapNbasesMin 10 \
  --alignSplicedMateMapLminOverLmate 0.5 \
  --alignSJstitchMismatchNmax 5 -1 5 5 \
  --chimSegmentMin 10 \
  --chimOutType WithinBAM HardClip \
  --chimJunctionOverhangMin 10 \
  --chimScoreDropMax 30 \
  --chimScoreJunctionNonGTAG 0 \
  --chimScoreSeparation 1 \
  --chimSegmentReadGapMax 3 \
  --chimMultimapNmax 50
```

---

## Step 4B. Run STAR alignment for STAR-Fusion-compatible chimeric detection

For STAR-Fusion-compatible alignment, run STAR with the chimeric-junction output settings expected by STAR-Fusion.

```bash
STAR \
  --genomeDir ~/references/GRCh38/gencode_v46/star \
  --readFilesIn ~/GSE112705/fastq/SRRXXXXXXX.fastq.gz \
  --readFilesCommand zcat \
  --runThreadN 8 \
  --outFileNamePrefix ~/GSE112705/star_prep/SRRXXXXXXX/SRRXXXXXXX. \
  --twopassMode Basic \
  --outReadsUnmapped None \
  --outSAMstrandField intronMotif \
  --outSAMunmapped Within \
  --chimSegmentMin 12 \
  --chimJunctionOverhangMin 8 \
  --chimOutJunctionFormat 1 \
  --alignSJDBoverhangMin 10 \
  --alignMatesGapMax 100000 \
  --alignIntronMax 100000 \
  --alignSJstitchMismatchNmax 5 -1 5 5 \
  --chimMultimapScoreRange 3 \
  --chimScoreJunctionNonGTAG -4 \
  --chimMultimapNmax 20 \
  --chimNonchimScoreDropMin 10 \
  --peOverlapNbasesMin 12 \
  --peOverlapMMp 0.1 \
  --alignInsertionFlush Right \
  --alignSplicedMateMapLminOverLmate 0 \
  --alignSplicedMateMapLmin 30 \
  --chimOutType Junctions \
  --quantMode GeneCounts \
  --outSAMtype BAM SortedByCoordinate
```

This configuration produces sorted BAM files, splice-junction files, gene-count summaries, and chimeric-junction outputs required for STAR-Fusion-compatible downstream analysis. The `--twopassMode Basic` option enables two-pass alignment, allowing STAR to detect splice junctions during an initial alignment pass and incorporate these junctions into a second alignment pass to improve mapping accuracy. `--outReadsUnmapped None` prevents unmapped reads from being written as separate output files, while `--readFilesCommand zcat` allows compressed FASTQ files to be processed directly without decompression. Strand information is retained using `--outSAMstrandField intronMotif`, and unmapped reads are preserved within the BAM file using `--outSAMunmapped Within`. Fusion detection sensitivity is controlled through `--chimSegmentMin 12`, which requires a minimum of 12 aligned bases per chimeric segment, and `--chimJunctionOverhangMin 8`, which requires at least 8 bases flanking a candidate chimeric breakpoint. Chimeric junction output is generated using `--chimOutJunctionFormat 1`, while `--alignSJDBoverhangMin 10` specifies the minimum overhang length required for splice-junction database alignment. Long-range fusion events are accommodated through `--alignMatesGapMax 100000` and `--alignIntronMax 100000`, allowing large genomic distances between aligned segments. Junction stitching behavior is controlled by `--alignSJstitchMismatchNmax 5 -1 5 5`, whereas `--chimMultimapScoreRange 3`, `--chimScoreJunctionNonGTAG -4`, `--chimMultimapNmax 20`, and `--chimNonchimScoreDropMin 10` regulate multimapping behavior, non-canonical splice junction penalties, and discrimination between chimeric and non-chimeric alignments. Parameters `--peOverlapNbasesMin 12` and `--peOverlapMMp 0.1` define overlap and mismatch tolerances for paired-end reads, while `--alignInsertionFlush Right` controls insertion placement relative to splice junctions. Minimum alignment requirements for spliced reads are specified through `--alignSplicedMateMapLminOverLmate 0` and `--alignSplicedMateMapLmin 30`. Chimeric junctions are exported using `--chimOutType Junctions`, gene-level counts are generated using `--quantMode GeneCounts`, and coordinate-sorted BAM files are produced using `--outSAMtype BAM SortedByCoordinate`, making all outputs directly compatible with STAR-Fusion and nf-core/rnafusion downstream workflows.

---

## Step 4C. Index coordinate-sorted BAM files

Index the coordinate-sorted BAM files using Samtools.

```bash
samtools index ~/GSE112705/star_prep/SRRXXXXXXX_Aligned.sortedByCoord.out.bam
```

The resulting `.bai` files are required for efficient genomic access and for compatibility with downstream workflow steps.

---

## Step 5. Prepare the nf-core/rnafusion samplesheet

Create a CSV file with the following structure:

```csv
sample,strandedness,bam,bai,junctions,splice_junctions
SRR6939927,forward,/path/to/SRR6939927.bam,/path/to/SRR6939927.bai,/path/to/Chimeric.out.junction,/path/to/SJ.out.tab
```

---

## Step 6. Build nf-core/rnafusion references

```bash
nextflow run nf-core/rnafusion \
  --build_references --all \
  --genomes_base ~/references/rnafusion \
  --outdir ~/references/rnafusion
```

This step prepares the required reference resources for the selected fusion callers. Because reference construction is computationally intensive, it should be completed before sample-level analysis and reused across runs when the same genome build is used.

---

## Step 7. Run nf-core/rnafusion

The pipeline is run using the samplesheet generated in the previous step.

```bash
nextflow run nf-core/rnafusion -r 4.1.0 \
  -profile singularity \
  --tools arriba,starfusion \
  --skip_qc \
  --input samplesheet_bam.csv \
  --genomes_base references \
  --outdir results \
  -work-dir work \
  -resume
```

The `--tools` argument specifies the fusion callers used in the analysis; here, Arriba and STAR-Fusion are selected. The `--genomes_base` argument defines the reference directory, `--outdir` defines the output directory, and `-resume` allows interrupted or repeated runs to continue from completed workflow steps.

---
## 3.2 Fusion Partner Identification
This section identifies recurrent fusion partner genes from Arriba fusion-calling results and generates a tumor-exclusive fusion partner gene set for downstream translatome analyses.

The complete implementation is provided in:

```r
section3.2_partner_extraction.R
```

All example input and output files used in this workflow are available in the repository. Users may replace these files with their own Arriba outputs while preserving the same structure.

---

## Step 8. Configure the R Environment and Import Arriba Outputs

The workflow begins by loading the required R packages and defining the directory containing the Arriba fusion-calling results. Each Arriba result file is associated with a sequencing run identifier, patient identifier, and tissue type (tumor or normal).

```r
library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(ggplot2)

arriba_dir <- "path/to/arriba_outputs"
```

A lookup table is then used to map sequencing runs to patient and tissue information, allowing downstream comparisons between tumor and matched normal samples.

---

## Step 9. Merge Fusion Events Across Samples

All Arriba fusion-calling results are imported and merged into a unified event-level dataset.

```r
read_arriba_events <- function(arriba_dir,
                               pattern = "\\.fusions\\.discarded\\.tsv$") {

  tsv_files <- list.files(arriba_dir,
                          pattern = pattern,
                          full.names = TRUE)

  event_list <- lapply(tsv_files, function(tsv) {

    df <- read_tsv(
      tsv,
      comment = "",
      show_col_types = FALSE,
      col_types = cols(.default = "c")
    )

    colnames(df)[1] <- sub("^#", "", colnames(df)[1])

    sample_id <- sub(pattern, "", basename(tsv))

    df %>% mutate(sample_id = sample_id, .before = 1)

  })

  bind_rows(event_list)
}

arriba_events <- read_arriba_events(arriba_dir)
```

All columns are initially imported as character variables to avoid type conflicts caused by Arriba missing-value symbols (`.`). Read-support metrics are subsequently converted to numeric values and annotated with patient and tissue information.

The resulting dataset contains one row per fusion event across the entire cohort.

---

## Step 10. Generate Fusion Partner Gene Summaries

Each fusion event contains a 5′ fusion partner (`gene1`) and a 3′ fusion partner (`gene2`). These partners are converted into a long-format representation and summarized at the gene level.

```r
partner_summary <- partner_long %>%
  group_by(gene_symbol) %>%
  summarise(
    n_samples        = n_distinct(sample_id),
    n_events         = n(),
    n_samples_tumor  = n_distinct(sample_id[tissue == "tumor"]),
    n_samples_normal = n_distinct(sample_id[tissue == "normal"]),
    .groups = "drop"
  )
```

For each fusion partner gene, the workflow calculates recurrence statistics, confidence-level distributions, reading-frame information, breakpoint characteristics, and read-support metrics.

---

## Step 11. Identify Tumor-Exclusive Fusion Partners

Fusion partner genes are filtered to retain only genes detected in tumor samples and absent from matched normal samples.

```r
partner_tumor_exclusive <- partner_summary %>%
  filter(
    n_samples_tumor >= 1,
    n_samples_normal == 0
  )
```

An additional recurrence threshold is applied to prioritize highly recurrent tumor-specific fusion partners. The workflow also evaluates multiple recurrence cutoffs and generates summary plots to assist threshold selection.

### Output Files

| File                                 | Description                        |
| ------------------------------------ | ---------------------------------- |
| arriba_events_long.csv               | Unified event-level fusion dataset |
| fusion_partner_events_long.csv       | Long-format fusion partner table   |
| fusion_partner_summary_full.csv      | Complete fusion partner summary    |
| fusion_partner_tumor_exclusive.csv   | Tumor-exclusive fusion partners    |
| fusion_partner_summary_filtered.csv  | Recurrence-filtered partners       |
| tumor_exclusive_recurrence_sweep.csv | Recurrence threshold statistics    |
| tumor_exclusive_recurrence_sweep.png | Recurrence threshold visualization |

---

# 3.4 Differential Translation Efficiency Analysis

This section performs differential translation efficiency analysis using matched RNA-seq and Ribo-seq datasets through the anota2seq framework.

The complete implementation is provided in:

```r
section3.4_anota2seq_run.R
```

All example input and output files are available within the repository.

---

## Step 13. Generate Ribo-seq and RNA-seq Quantification Data

Gene-level RNA-seq and Ribo-seq count matrices are generated following alignment and quantification. These count matrices serve as the input for differential translation efficiency analysis.

---

## Step 14. Prepare Count Matrices for anota2seq

The combined count matrix is imported into R and subjected to quality-control filtering.

```r
dat <- read.csv(
  "GSE112705_RPF_RNA_readCounts.csv",
  header = TRUE,
  stringsAsFactors = FALSE
)

dat$Gene_name <- make.unique(as.character(dat$Gene_name))
rownames(dat) <- dat$Gene_name
```

Genes containing zero counts across samples are removed prior to analysis.

The matrix is then separated into translated mRNA (RPF) and total mRNA (RNA) datasets.

```r
my_data_P <- dat_sorted[, grep("RPF", colnames(dat_sorted))]
my_data_T <- dat_sorted[, grep("RNA", colnames(dat_sorted))]
```

Phenotype labels and batch information are automatically extracted from sample identifiers.

---

## Step 15. Construct the anota2seq Dataset

The RNA-seq and Ribo-seq matrices are combined into a single anota2seq object.

```r
ads <- anota2seqDataSetFromMatrix(
  dataP = as.matrix(my_data_P),
  dataT = as.matrix(my_data_T),
  phenoVec = myPheno,
  batchVec = myBatch,
  dataType = "RNAseq",
  normalize = TRUE,
  transformation = "TMM-log2"
)
```

Normalization is performed using the recommended TMM-log2 strategy.

---

## Step 16. Perform Differential Translation Efficiency Analysis

Differential translational regulation is evaluated using the Analysis of Partial Variance (APV) framework.

```r
ads <- anota2seqAnalyze(
  Anota2seqDataSet = ads,
  contrasts = myContrast,
  analysis = c(
    "translation",
    "buffering",
    "translated mRNA",
    "total mRNA"
  )
)
```

The workflow simultaneously evaluates:

* Translation efficiency
* Buffering effects
* Translated mRNA abundance
* Total mRNA abundance

---

## Step 17. Classify Genes into Regulatory Modes

Genes are assigned to regulatory categories according to the anota2seq hierarchical framework.

```r
runAds <- anota2seqRegModes(runAds)

dataOut <- anota2seqGetOutput(
  runAds,
  output = "singleDf",
  selContrast = 1
)
```

Genes are classified into:

* Translation
* Buffering
* Abundance
* Background

### Output Files

| File                                    | Description                            |
| --------------------------------------- | -------------------------------------- |
| GSE112705_anota2seq_Results_Cleaned.csv | Final differential translation results |
| GSE112705_PreAnalysis_PCA.pdf           | PCA quality-control report             |
| GSE112705_RegulatoryModes_Scatter.png   | Regulatory mode visualization          |

---

# 3.5 Translatome-Guided Prioritization of Fusion Partner Genes

This section integrates tumor-exclusive fusion partner genes with translation-efficiency results to identify fusion-associated genes exhibiting evidence of translational regulation.

The complete implementation is provided in:

```r
section3.5_partner_translation.R
```

All example input and output files are available within the repository.

---

## Step 18. Integrate Fusion Partners with Translation Efficiency Results

Tumor-exclusive fusion partners generated in Section 3.2 are merged with the translation-efficiency results generated in Section 3.4.

```r
partner_TE <- partner_raw %>%
  inner_join(anota_raw, by = "gene_symbol")
```

The resulting table combines fusion recurrence statistics with translational regulation measurements.

A fusion-partner-only translation-efficiency table is then generated.

```r
anota_partners_only <- partner_TE %>%
  select(
    gene_symbol,
    translatedmRNA.apvEff,
    totalmRNA.apvEff,
    translation.apvEff,
    buffering.apvEff,
    singleRegMode
  )
```

Regulatory-mode distributions are subsequently summarized.

```r
mode_summary <- partner_TE %>%
  count(recurrence_tier, singleRegMode)
```

---

## Step 19. Visualize Translational Regulation of Fusion Partner Genes

Fusion partner genes are visualized using an anota2seq-style fold-change scatter plot.

```r
fig_A <- ggplot(
  fp_df,
  aes(deltaT, deltaP, colour = group)
) +
  geom_point() +
  theme_classic()
```

The horizontal axis represents total mRNA abundance changes, whereas the vertical axis represents translated mRNA abundance changes. Genes are coloured according to their assigned regulatory mode and direction of regulation.

This visualization facilitates the identification of:

* Translational activation
* Translational repression
* Buffering effects
* Abundance-driven regulation

### Output Files

| File                         | Description                                                  |
| ---------------------------- | ------------------------------------------------------------ |
| fusion_partner_TE_table.csv  | Integrated fusion partner and translation-efficiency results |
| anota2seq_partners_only.csv  | Translation-efficiency results restricted to fusion partners |
| partner_mode_summary.csv     | Regulatory mode summary                                      |
| fig_1_partner_FC_scatter.pdf | Publication-quality figure                                   |
| fig_1_partner_FC_scatter.png | PNG version of the figure                                    |
E_table.csv`
* `anota2seq_partners_only.csv`
* `partner_mode_summary.csv`
* `fig_1_partner_FC_scatter.pdf`
* `fig_1_partner_FC_scatter.png`
