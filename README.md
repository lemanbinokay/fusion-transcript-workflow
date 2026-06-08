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

### Step 8. Configure the R Environment and Import Arriba Outputs

The analysis begins by loading the required R packages and defining the directory containing the Arriba fusion-calling results. All Arriba output files should be stored within a single directory. A sample metadata table is then constructed to map sequencing run identifiers to patient IDs and tissue types, enabling downstream comparison of tumor and normal samples.

The complete implementation of this step is provided in:

`section3.2_partner_extraction.R`

### Step 9. Merge Fusion Events Across Samples

All Arriba fusion call files are imported and merged into a unified event-level dataset. Fusion-support metrics, breakpoint annotations, confidence levels, and additional Arriba annotations are standardized across samples. Each fusion event is subsequently linked to the corresponding patient and tissue information through the sample metadata table.

The resulting dataset contains one row per fusion event and serves as the foundation for all downstream fusion partner analyses.

### Step 10. Generate Fusion Partner Gene Summaries

Fusion partner genes are extracted from both the 5′ and 3′ breakpoints of each fusion event and transformed into a long-format partner table. Information from all events involving the same gene is then aggregated to generate gene-level summaries.

For each fusion partner gene, recurrence statistics, confidence metrics, breakpoint characteristics, read-support information, and tissue distribution are calculated. This step provides a comprehensive overview of fusion partner occurrence across the entire cohort.

### Step 11. Identify Tumor-Exclusive Fusion Partners

Fusion partner genes are filtered to retain only those detected in tumor samples and absent from all matched normal samples. These tumor-exclusive partners constitute the candidate gene set used for downstream translational analyses.

An optional recurrence filter can be applied to retain only genes observed in multiple tumor samples. To facilitate threshold selection, recurrence summaries and visualization plots are generated across a range of recurrence cutoffs.

**Output files**

* `arriba_events_long.csv`
* `fusion_partner_events_long.csv`
* `fusion_partner_summary_full.csv`
* `fusion_partner_tumor_exclusive.csv`
* `fusion_partner_summary_filtered.csv`
* `tumor_exclusive_recurrence_sweep.csv`
* `tumor_exclusive_recurrence_sweep.png`

---

## 3.4 Differential Translation Efficiency Analysis

### Step 13. Generate Ribo-seq and RNA-seq Quantification Data

Translation efficiency analysis requires matched Ribo-seq and RNA-seq datasets processed through a standardized quantification workflow. Libraries are organized using a sample sheet that specifies sample identifiers, sequencing files, strandedness information, and assay type. Following alignment and quantification, gene-level count matrices are generated for both translated mRNA and total mRNA fractions.

### Step 14. Prepare Count Matrices for anota2seq

The combined count matrix is imported into R and subjected to quality-control filtering. Gene identifiers are standardized, non-numeric values are removed, and genes containing zero counts across samples are excluded. The filtered matrix is subsequently separated into translated-mRNA (RPF) and total-mRNA (RNA) components.

Sample conditions and batch information are extracted directly from the library identifiers to support downstream statistical modelling.

The complete implementation of this step is provided in:

`section3.4_anota2seq_run.R`

### Step 15. Construct the anota2seq Dataset

The translated-mRNA matrix, total-mRNA matrix, phenotype labels, and batch information are combined into a single anota2seq object. Normalization is performed using the recommended TMM-log2 transformation, providing a unified framework for downstream translational regulation analyses.

### Step 16. Perform Differential Translation Efficiency Analysis

Differential translational regulation is assessed using the Analysis of Partial Variance (APV) framework implemented in anota2seq. Translation, buffering, translated-mRNA abundance, and total-mRNA abundance analyses are performed simultaneously, allowing comprehensive characterization of gene regulation between tumor and normal tissues.

### Step 17. Classify Genes into Regulatory Modes

Genes are assigned to regulatory categories according to the anota2seq hierarchical classification framework. Each gene is classified into one of four regulatory modes:

* Translation
* Buffering
* Abundance
* Background

The resulting dataset provides the basis for identifying translationally regulated fusion partner genes.

**Output files**

* `GSE112705_anota2seq_Results_Cleaned.csv`
* `GSE112705_PreAnalysis_PCA.pdf`
* `GSE112705_RegulatoryModes_Scatter.png`

---

## 3.5 Translatome-Guided Prioritization of Fusion Partner Genes

### Step 18. Integrate Fusion Partners with Translation Efficiency Results

Tumor-exclusive fusion partner genes are integrated with the anota2seq results to identify fusion-associated genes displaying evidence of translational regulation. Gene-level translation efficiency statistics, abundance changes, buffering effects, and regulatory mode assignments are merged into a unified dataset.

Summary tables are generated to facilitate downstream interpretation of translationally regulated fusion partners and to quantify the distribution of regulatory modes across recurrence categories.

The complete implementation of this step is provided in:

`section3.5_partner_translation.R`

### Step 19. Visualize Translational Regulation of Fusion Partner Genes

Fusion partner genes are visualized using an anota2seq-style fold-change scatter plot. The horizontal axis represents changes in total mRNA abundance, whereas the vertical axis represents changes in translated mRNA abundance. Genes are coloured according to their assigned regulatory mode and direction of regulation.

<img width="1920" height="1080" alt="CHAPTER_FIGURE" src="https://github.com/user-attachments/assets/8821197d-497c-4365-9a0b-2dde0d6a8756" />


This visualization enables rapid identification of fusion partner genes exhibiting translational activation, translational repression, buffering effects, or abundance-driven regulation.

**Output files**

* `fusion_partner_TE_table.csv`
* `anota2seq_partners_only.csv`
* `partner_mode_summary.csv`
* `fig_1_partner_FC_scatter.pdf`
* `fig_1_partner_FC_scatter.png`
