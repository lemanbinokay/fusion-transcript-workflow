# ---- 1. Setup ---------------------------------------------------------------

required_pkgs <- c("dplyr", "readr", "tidyr", "stringr", "ggplot2",
                   "ggrepel", "patchwork", "scales", "broom")
missing_pkgs  <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(missing_pkgs) > 0) install.packages(missing_pkgs)

library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(scales)

arriba_dir <- "C:/Users/baris/Desktop/arriba"

# Create and set a distinct working directory for these results
work_dir <- file.path(arriba_dir, "tumor_exclusive_analysis")
if (!dir.exists(work_dir)) {
  dir.create(work_dir, recursive = TRUE)
  message("Created new output directory: ", work_dir)
}
setwd(work_dir)

# Strictly enforce loading of the cleaned anota2seq data to guarantee pipeline consistency
anota_path <- file.path(arriba_dir, "GSE112705_anota2seq_Results_Cleaned.csv")
stopifnot("Cleaned anota2seq file missing" = file.exists(anota_path))
message("Reading cleaned anota2seq output from: ", anota_path)

# Target the tumor-exclusive file without recurrence filtering
partner_path <- file.path(arriba_dir, "fusion_partner_tumor_exclusive.csv")
stopifnot("Tumor-exclusive partner file missing" = file.exists(partner_path))


# ---- 2. Load and harmonise the two tables -----------------------------------

anota_raw   <- read_csv(anota_path,   show_col_types = FALSE)
partner_raw <- read_csv(partner_path, show_col_types = FALSE)

# RE-INJECT METADATA: Restore the categorical recurrence_tier column to 
# maintain compatibility with downstream statistical grouping.
if (!"recurrence_tier" %in% colnames(partner_raw)) {
  partner_raw <- partner_raw %>% mutate(recurrence_tier = "tumor_exclusive")
}

if ("Gene" %in% colnames(anota_raw) && !"gene_symbol" %in% colnames(anota_raw)) {
  anota_raw <- anota_raw %>% rename(gene_symbol = Gene)
}

# Ensure genes that did not pass filters are properly classified as background
anota_raw <- anota_raw %>%
  mutate(singleRegMode = ifelse(is.na(singleRegMode) | singleRegMode == "",
                                "background",
                                str_trim(singleRegMode)))

need_cols <- c("gene_symbol",
               "translatedmRNA.apvEff", "translatedmRNA.apvRvmPAdj",
               "totalmRNA.apvEff",      "totalmRNA.apvRvmPAdj",
               "translation.apvEff",    "translation.apvRvmPAdj",
               "buffering.apvEff",      "buffering.apvRvmPAdj",
               "singleRegMode")
missing <- setdiff(need_cols, colnames(anota_raw))
if (length(missing) > 0) {
  stop("anota2seq output is missing columns: ", paste(missing, collapse = ", "))
}


# ---- 3. Inner join: fusion partners x anota2seq -----------------------------

partner_TE <- partner_raw %>%
  inner_join(anota_raw, by = "gene_symbol")

write_csv(partner_TE, file.path(work_dir, "fusion_partner_TE_table.csv"))

anota_partners_only <- partner_TE %>%
  select(gene_symbol,
         translatedmRNA.apvEff, translatedmRNA.apvRvmPAdj,
         totalmRNA.apvEff,      totalmRNA.apvRvmPAdj,
         translation.apvEff,    translation.apvRvmPAdj,
         buffering.apvEff,      buffering.apvRvmPAdj,
         singleRegMode)

write_csv(anota_partners_only, file.path(work_dir, "anota2seq_partners_only.csv"))

message("Fusion partners with anota2seq results: ", nrow(partner_TE),
        " of ", nrow(partner_raw),
        " filtered partners (",
        round(100 * nrow(partner_TE) / nrow(partner_raw), 1), "%)")
message("anota2seq_partners_only.csv written: ",
        nrow(anota_partners_only), " genes, ",
        ncol(anota_partners_only), " columns")


# ---- 4. Per-tier mode summary (text output) ---------------------------------

mode_summary <- partner_TE %>%
  count(recurrence_tier, singleRegMode) %>%
  pivot_wider(names_from = singleRegMode, values_from = n, values_fill = 0) %>%
  arrange(recurrence_tier)

write_csv(mode_summary, file.path(work_dir, "partner_mode_summary.csv"))
message("Per-tier mode summary:")
print(mode_summary)


# ---- 5. Shared plotting helpers --------------------------------------------

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

make_plot_df <- function(df) {
  df %>%
    mutate(deltaP = translatedmRNA.apvEff,
           deltaT = totalmRNA.apvEff,
           group  = direction_label(singleRegMode, deltaP, deltaT))
}

bg_df <- make_plot_df(anota_raw)
fp_df <- make_plot_df(partner_TE)

ax_lim <- 5


# ---- 6. Figure 1: partner-only fold-change scatter --------------------------


fp_counts <- fp_df %>% 
  count(group, .drop = FALSE) %>%
  mutate(label = sprintf("%s (n=%d)", group, n))
fp_legend_labels <- setNames(fp_counts$label, fp_counts$group)

fig_A <- ggplot(fp_df, aes(deltaT, deltaP, colour = group)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey60") +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey60") +
  geom_point(data = filter(fp_df, group == "Background"),
             colour = "#CCCCCC", size = 1.0, alpha = 0.55) +
  geom_point(data = filter(fp_df, group != "Background"),
             size = 1.4, alpha = 0.85) +
  scale_colour_manual(values = mode_colors,
                      labels = fp_legend_labels,
                      breaks = mode_levels,
                      drop   = FALSE) +
  coord_cartesian(xlim = c(-ax_lim, ax_lim), ylim = c(-ax_lim, ax_lim)) +
  labs(x = "Total mRNA (Log2FC)",
       y = "Translated mRNA (Log2FC)",
       colour = NULL,
       title  = "Tumor vs Normal: fusion-partner translatome") +
  theme_classic(base_size = 12) +
  theme(legend.position = c(0.02, 0.98),
        legend.justification = c(0, 1),
        legend.background = element_rect(fill = scales::alpha("white", 0.8),
                                         colour = NA),
        legend.text = element_text(size = 9),
        plot.title = element_text(face = "bold", hjust = 0.5))

ggsave(file.path(work_dir, "fig_1_partner_FC_scatter.pdf"),
       fig_A, width = 7, height = 7)
ggsave(file.path(work_dir, "fig_1_partner_FC_scatter.png"),
       fig_A, width = 7, height = 7, dpi = 300)