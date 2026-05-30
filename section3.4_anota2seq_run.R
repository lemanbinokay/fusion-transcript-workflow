
################################################################################
# Comprehensive anota2seq Pipeline: GSE112705 (Tumor vs Normal)
################################################################################

library(anota2seq)
library(limma)
library(edgeR)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(tools)
library(stringr)

# ---------------------------------------------------------
# 1. Directory Setup & Data Ingestion
# ---------------------------------------------------------
# Force R to your specific working directory
base_dir <- "C:/Users/baris/Desktop/arriba"
setwd(base_dir)
message("Base directory set to: ", getwd())

# Explicitly check if the file exists before proceeding
target_file <- "GSE112705_RPF_RNA_readCounts.csv"
if(!file.exists(target_file)) {
  stop("CRITICAL ERROR: Cannot find ", target_file, " in ", base_dir, ". Check the spelling or location.")
}

# Create a clean output folder inside the arriba directory and move into it
out_dir <- "GSE112705_anota2seq_Comprehensive_Output"
if(!dir.exists(out_dir)) dir.create(out_dir)
setwd(out_dir)
message("Output directory set to: ", getwd())

# Read the file (stepping back one directory level to the base folder)
dat <- read.csv(paste0("../", target_file), header = TRUE, stringsAsFactors = FALSE)

# ---------------------------------------------------------
# 2. Strict Cleaning & Matrix Alignment
# ---------------------------------------------------------
dat$Gene_name <- make.unique(as.character(dat$Gene_name))
rownames(dat) <- dat$Gene_name


# Isolate counts and force numeric conversion
dat_counts <- dat[, !(names(dat) %in% c("geneID", "Gene_name"))]
dat_counts <- as.data.frame(lapply(dat_counts, function(x) as.numeric(as.character(x))))
rownames(dat_counts) <- rownames(dat)
dat_counts <- na.omit(dat_counts)

# Strict Zero Filtering for RVM assumptions
non_zero_rows <- apply(dat_counts, 1, function(x) all(x > 0))
dat_filtered <- dat_counts[non_zero_rows, ]

if(nrow(dat_filtered) == 0) stop("CRITICAL: Filtering removed all genes.")

dat_sorted <- dat_filtered[, order(colnames(dat_filtered))]
my_data_P <- dat_sorted[, grep("RPF", colnames(dat_sorted))]
my_data_T <- dat_sorted[, grep("RNA", colnames(dat_sorted))]

if(!all(gsub("RPF", "", colnames(my_data_P)) == gsub("RNA", "", colnames(my_data_T)))) {
  stop("CRITICAL: RPF and RNA sample matrices are misaligned.")
}

get_condition <- function(x) {
  if (grepl("normal", x, ignore.case = TRUE)) return("Normal")
  if (grepl("tumor", x, ignore.case = TRUE)) return("Tumor")
}
myPheno <- sapply(colnames(my_data_P), get_condition)
myBatch <- sapply(strsplit(colnames(my_data_P), "[\\.-]"), `[`, 1)

# ---------------------------------------------------------
# 3. Pre-Analysis PCA (Manual Variance Control)
# ---------------------------------------------------------
pca_norm <- voom(calcNormFactors(DGEList(as.matrix(dat_filtered))))$E
sd_vals <- apply(pca_norm, 1, sd)
pca_filtered <- pca_norm[sd_vals > quantile(sd_vals, 0.25), ]

pca_out <- prcomp(t(pca_filtered))
anot_pca <- data.frame(
  sample = colnames(dat_filtered),
  condition = sapply(colnames(dat_filtered), get_condition),
  modality = ifelse(grepl("RPF", colnames(dat_filtered)), "RPF", "RNA"),
  batch = sapply(strsplit(colnames(dat_filtered), "[\\.-]"), `[`, 1)
)
rownames(anot_pca) <- anot_pca$sample
pca_plot_df <- merge(pca_out$x, anot_pca, by = "row.names")

pdf("GSE112705_PreAnalysis_PCA.pdf", width = 12, height = 6)
grid.arrange(
  ggplot(pca_plot_df, aes(x = PC1, y = PC2, shape = modality, col = condition)) +
    geom_point(size = 3) + theme_bw() + ggtitle("PC1 vs PC2 Variance"),
  ggplot(pca_plot_df, aes(x = PC1, y = PC3, shape = modality, col = condition)) +
    geom_point(size = 3) + theme_bw() + ggtitle("PC1 vs PC3 Variance"),
  ncol = 2
)
dev.off()

# ---------------------------------------------------------
# 4. anota2seq Initialization & Execution
# ---------------------------------------------------------
ads <- anota2seqDataSetFromMatrix(
  dataP           = as.matrix(my_data_P),
  dataT           = as.matrix(my_data_T),
  phenoVec        = myPheno,
  batchVec        = myBatch,
  dataType        = "RNAseq",
  normalize       = TRUE,
  transformation  = "TMM-log2",
  filterZeroGenes = TRUE,
  varCutOff       = NULL
)


phenoLev <- levels(as.factor(myPheno))
myContrast <- matrix(nrow = length(phenoLev), ncol = 1, dimnames = list(phenoLev, "Tumor_vs_Normal"))
myContrast[which(phenoLev == "Normal"), 1] <- -1
myContrast[which(phenoLev == "Tumor"), 1] <- 1

# Explicit step-wise execution to ensure diagnostic outputs

ads <- anota2seqAnalyze(
  Anota2seqDataSet = ads,
  contrasts        = myContrast,
  analysis = c("translation", "buffering",
               "translated mRNA", "total mRNA")
)

runAds <- anota2seqRun(
  Anota2seqDataSet = ads, 
  contrasts = myContrast, 
  performQC = TRUE,
  onlyGroup = FALSE, 
  performROT = TRUE, 
  generateSingleGenePlots = TRUE, # Set to TRUE if you want per-gene plots
  analyzeBuffering = TRUE, 
  analyzemRNA = TRUE,
  useRVM = TRUE, 
  correctionMethod = "BH", 
  useProgBar = TRUE,
  thresholds = list(
    maxPAdj = 0.25,
    deltaP = log2(1.2), 
    deltaT = log2(1.2),
    deltaPT = log2(1.2),
    deltaTP = log2(1.2), 
    maxSlopeTranslation = 1.5,
    minSlopeTranslation = -0.5,
    minSlopeBuffering = -1.5, 
    maxSlopeBuffering = 0.5
  )
)


# ---------------------------------------------------------
# 5. Generate QC plots and fold-change plots
# ---------------------------------------------------------

anota2seqPlotPvalues(runAds, selContrast = 1, plotToFile = FALSE,
                     contrastName = "Tumor vs Normal")
anota2seqPlotFC(     runAds, selContrast = 1, plotToFile = FALSE,
                     contrastName = "Tumor vs Normal")


# ---------------------------------------------------------
# 6. Run regulatory mode classification
# ---------------------------------------------------------
runAds <- anota2seqRegModes(runAds)


# ---------------------------------------------------------
# 7. Extract final results
# ---------------------------------------------------------
dataOut <- anota2seqGetOutput(runAds, output = "singleDf", selContrast = 1)

# Save results to CSV
write.csv(dataOut,
          file = "GSE112705_anota2seq_Results_Cleaned.csv",
          row.names = FALSE)

# ---------------------------------------------------------
# 8. Dynamic Regulatory Mode Scatter Plot
# ---------------------------------------------------------

library(scales)

df <- dataOut

# Define ordered mode levels and colors
mode_levels <- c("Translation up", "Translation down",
                 "Buffering up",   "Buffering down",
                 "Abundance up",   "Abundance down",
                 "Background")

mode_colors <- c(
  "Translation up"   = "#E8765B",
  "Translation down" = "#A6231C",
  "Buffering up"     = "#9DC4E0",
  "Buffering down"   = "#1F4E7A",
  "Abundance up"     = "#A6D49F",
  "Abundance down"   = "#2F7A36",
  "Background"       = "#BDBDBD"
)

# Assign directional group based on regulatory mode
direction_label <- function(mode, deltaP, deltaT) {
  out <- rep("Background", length(mode))
  out[mode == "translation" & deltaP >= 0] <- "Translation up"
  out[mode == "translation" & deltaP <  0] <- "Translation down"
  out[mode == "abundance"   & deltaP >= 0] <- "Abundance up"
  out[mode == "abundance"   & deltaP <  0] <- "Abundance down"
  out[mode == "buffering"   & deltaT >= 0] <- "Buffering up"
  out[mode == "buffering"   & deltaT <  0] <- "Buffering down"
  factor(out, levels = mode_levels)
}

df <- df %>%
  mutate(singleRegMode = ifelse(is.na(singleRegMode) | singleRegMode == "",
                                "background", str_trim(singleRegMode)),
         deltaP = translatedmRNA.apvEff,
         deltaT = totalmRNA.apvEff,
         group  = direction_label(singleRegMode, deltaP, deltaT))

# Build legend labels with counts
grp_counts <- df %>%
  count(group, .drop = FALSE) %>%
  mutate(label = sprintf("%s (n=%d)", group, n))
legend_labels <- setNames(grp_counts$label, grp_counts$group)

ax_lim <- 5

p <- ggplot(df, aes(deltaT, deltaP, colour = group)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey60") +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey60") +
  geom_point(data = filter(df, group == "Background"),
             colour = "#CCCCCC", size = 1.0, alpha = 0.55) +
  geom_point(data = filter(df, group != "Background"),
             size = 1.4, alpha = 0.85) +
  scale_colour_manual(values = mode_colors,
                      labels = legend_labels,
                      breaks = mode_levels,
                      drop   = FALSE) +
  coord_cartesian(xlim = c(-ax_lim, ax_lim), ylim = c(-ax_lim, ax_lim)) +
  labs(x = "Total mRNA (Log2FC)",
       y = "Translated mRNA (Log2FC)",
       colour = NULL,
       title  = "GSE112705 Tumor vs Normal: Translatome Analysis") +
  theme_classic(base_size = 12) +
  theme(legend.position = c(0.02, 0.98),
        legend.justification = c(0, 1),
        legend.background = element_rect(fill = scales::alpha("white", 0.8),
                                         colour = NA),
        legend.text = element_text(size = 9),
        plot.title = element_text(face = "bold", hjust = 0.5))

ggsave("GSE112705_RegulatoryModes_Scatter.png",
       plot = p, width = 7, height = 7, dpi = 300)



## automatic graph throuh anota2seq
anota2seqPlotPvalues(runAds, selContrast = 1, plotToFile = FALSE,
                     contrastName = "Tumor vs Normal")
anota2seqPlotFC(     runAds, selContrast = 1, plotToFile = FALSE,
                     contrastName = "Tumor vs Normal")


