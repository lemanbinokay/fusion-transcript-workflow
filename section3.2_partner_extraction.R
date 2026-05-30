

# ---- 1. Setup ---------------------------------------------------------------

required_pkgs <- c("dplyr", "readr", "tidyr", "stringr", "ggplot2")
missing_pkgs  <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(missing_pkgs) > 0) install.packages(missing_pkgs)

library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(ggplot2)

arriba_dir <- "C:/Users/baris/Desktop/arriba"

if (!dir.exists(arriba_dir)) {
  stop("Arriba directory not found: ", arriba_dir)
}


# ---- 2. SRR -> tissue / patient lookup --------------------------------------

sample_lookup <- tribble(
  ~srr_id,        ~patient,  ~tissue,
  "SRR6939925",   "LC001",   "normal",
  "SRR6939927",   "LC001",   "tumor",
  "SRR6939929",   "LC033",   "normal",
  "SRR6939931",   "LC033",   "tumor",
  "SRR6939933",   "LC034",   "normal",
  "SRR6939935",   "LC034",   "tumor",
  "SRR6939937",   "LC501",   "normal",
  "SRR6939939",   "LC501",   "tumor",
  "SRR6939941",   "LC502",   "normal",
  "SRR6939943",   "LC502",   "tumor",
  "SRR6939945",   "LC505",   "normal",
  "SRR6939947",   "LC505",   "tumor",
  "SRR6939949",   "LC506",   "normal",
  "SRR6939951",   "LC506",   "tumor",
  "SRR6939954",   "LC507",   "normal",
  "SRR6939956",   "LC507",   "tumor",
  "SRR6939958",   "LC508",   "normal",
  "SRR6939960",   "LC508",   "tumor",
  "SRR6939962",   "LC509",   "normal",
  "SRR6939964",   "LC509",   "tumor"
)

extract_srr <- function(x) str_extract(x, "SRR\\d{6,9}")


# ---- 3. Read all Arriba TSV files -------------------------------------------

read_arriba_events <- function(arriba_dir,
                               pattern = "\\.fusions\\.discarded\\.tsv$") {
  tsv_files <- list.files(arriba_dir, pattern = pattern, full.names = TRUE)
  if (length(tsv_files) == 0) {
    stop("No Arriba TSV files matched '", pattern, "' in: ", arriba_dir)
  }
  event_list <- lapply(tsv_files, function(tsv) {
    df <- read_tsv(tsv, comment = "", show_col_types = FALSE,
                   col_types = cols(.default = "c"))
    colnames(df)[1] <- sub("^#", "", colnames(df)[1])
    sample_id <- sub(pattern, "", basename(tsv))
    df %>% mutate(sample_id = sample_id, .before = 1)
  })
  bind_rows(event_list)
}

arriba_to_numeric <- function(x) {
  x <- ifelse(x %in% c(".", "", "NA"), NA_character_, x)
  suppressWarnings(as.numeric(x))
}

arriba_events <- read_arriba_events(arriba_dir)

arriba_events <- arriba_events %>%
  mutate(
    coverage1        = arriba_to_numeric(coverage1),
    coverage2        = arriba_to_numeric(coverage2),
    split_reads1     = arriba_to_numeric(split_reads1),
    split_reads2     = arriba_to_numeric(split_reads2),
    discordant_mates = arriba_to_numeric(discordant_mates)
  ) %>%
  mutate(srr_id = extract_srr(sample_id)) %>%
  left_join(sample_lookup, by = "srr_id")

unmatched <- arriba_events %>% filter(is.na(tissue)) %>% distinct(sample_id)
if (nrow(unmatched) > 0) {
  warning("Some sample_id values did not map to the SRR lookup: ",
          paste(unmatched$sample_id, collapse = ", "))
}

message("Total events read: ", nrow(arriba_events))
message("Samples represented: ", n_distinct(arriba_events$sample_id))
message("Tissue distribution (events per tissue):")
print(table(arriba_events$tissue, useNA = "ifany"))


# ---- 4. Resolve compound intergenic partner annotations ---------------------

resolve_partner <- function(gene_field) {
  primary <- str_split_fixed(gene_field, ",", 2)[, 1]
  primary <- str_replace(primary, "\\(.+\\)$", "")
  primary
}

arriba_events <- arriba_events %>%
  mutate(
    gene1_raw        = gene1,
    gene2_raw        = gene2,
    gene1            = resolve_partner(gene1),
    gene2            = resolve_partner(gene2),
    gene1_intergenic = str_detect(gene1_raw, ","),
    gene2_intergenic = str_detect(gene2_raw, ",")
  )


# ---- 5. Long-format partner table -------------------------------------------

partner_long <- bind_rows(
  arriba_events %>%
    transmute(
      sample_id, srr_id, patient, tissue,
      gene_symbol      = gene1,
      partner_position = "gene1",
      site             = site1,
      coverage         = coverage1,
      partner_other    = gene2,
      type, confidence, reading_frame,
      split_reads1, split_reads2, discordant_mates,
      intergenic       = gene1_intergenic
    ),
  arriba_events %>%
    transmute(
      sample_id, srr_id, patient, tissue,
      gene_symbol      = gene2,
      partner_position = "gene2",
      site             = site2,
      coverage         = coverage2,
      partner_other    = gene1,
      type, confidence, reading_frame,
      split_reads1, split_reads2, discordant_mates,
      intergenic       = gene2_intergenic
    )
)


# ---- 6. Per-partner summary (full cohort) -----------------------------------

safe_max <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) NA_real_ else max(x)
}

n_tumor_total  <- sample_lookup %>% filter(tissue == "tumor")  %>% nrow()
n_normal_total <- sample_lookup %>% filter(tissue == "normal") %>% nrow()
cohort_size    <- n_tumor_total + n_normal_total

# Recurrence cutoff: >= half of tumor samples
recurrence_cutoff <- ceiling(n_tumor_total / 2)   # = 5 for 10 tumor samples
message("Recurrence cutoff (>= half of ", n_tumor_total,
        " tumor samples): n_samples_tumor >= ", recurrence_cutoff)

partner_summary <- partner_long %>%
  group_by(gene_symbol) %>%
  summarise(
    n_samples              = n_distinct(sample_id),
    n_events               = n(),
    n_samples_tumor        = n_distinct(sample_id[tissue == "tumor"]),
    n_samples_normal       = n_distinct(sample_id[tissue == "normal"]),
    n_as_gene1             = sum(partner_position == "gene1"),
    n_as_gene2             = sum(partner_position == "gene2"),
    n_high_confidence      = sum(confidence == "high",   na.rm = TRUE),
    n_medium_confidence    = sum(confidence == "medium", na.rm = TRUE),
    n_low_confidence       = sum(confidence == "low",    na.rm = TRUE),
    n_inframe              = sum(reading_frame == "in-frame",     na.rm = TRUE),
    n_outofframe           = sum(reading_frame == "out-of-frame", na.rm = TRUE),
    n_frame_unknown        = sum(reading_frame == "." | is.na(reading_frame)),
    n_intergenic_flank     = sum(intergenic, na.rm = TRUE),
    sites_observed         = paste(sort(unique(site)), collapse = ","),
    types_observed         = paste(sort(unique(type)), collapse = ","),
    max_split_reads_total  = safe_max(split_reads1 + split_reads2),
    max_discordant_mates   = safe_max(discordant_mates),
    median_coverage        = median(coverage, na.rm = TRUE),
    samples                = paste(sort(unique(sample_id)), collapse = ","),
    .groups = "drop"
  ) %>%
  arrange(desc(n_samples), desc(n_events))


# ---- 7. Stage 1: Tumor-exclusive gate (n_normal == 0, n_tumor >= 1) ---------

partner_tumor_exclusive <- partner_summary %>%
  filter(
    n_samples_tumor  >= 1,
    n_samples_normal == 0
  ) %>%
  arrange(desc(n_samples_tumor), desc(n_events))

n_tumor_exclusive <- nrow(partner_tumor_exclusive)

message("\n--- Stage 1: Tumor-exclusive partners (n_normal == 0, n_tumor >= 1) ---")
message("Total: ", n_tumor_exclusive)
message("Tumor sample distribution:")
print(table(partner_tumor_exclusive$n_samples_tumor))


# ---- 8. Stage 2: Apply recurrence cutoff (n_samples_tumor >= 5) -------------

partner_summary_filtered <- partner_tumor_exclusive %>%
  filter(n_samples_tumor >= recurrence_cutoff) %>%
  arrange(desc(n_samples_tumor), desc(n_events))

n_filtered <- nrow(partner_summary_filtered)

message("\n--- Stage 2: After recurrence cutoff (n_samples_tumor >= ",
        recurrence_cutoff, ") ---")
message("Partners retained: ", n_filtered,
        " (", round(100 * n_filtered / n_tumor_exclusive, 1),
        "% of tumor-exclusive set)")
message("Partners dropped:  ", n_tumor_exclusive - n_filtered)


# ---- 9. Recurrence threshold sweep ------------------------------------------
# Show how many tumor-exclusive partners survive at every threshold n >= 1..10.
# The chosen cutoff (n=5) is highlighted.

sweep_thresholds <- seq_len(n_tumor_total)

recurrence_sweep <- tibble(
  min_tumor_samples   = sweep_thresholds,
  n_partners_retained = sapply(sweep_thresholds, function(t) {
    sum(partner_tumor_exclusive$n_samples_tumor >= t)
  }),
  pct_retained = sapply(sweep_thresholds, function(t) {
    100 * sum(partner_tumor_exclusive$n_samples_tumor >= t) /
      n_tumor_exclusive
  }),
  is_chosen = sweep_thresholds == recurrence_cutoff
)

message("\n--- Recurrence sweep (tumor-exclusive partners) ---")
print(recurrence_sweep, n = Inf)


# ---- 10. Sweep plot ----------------------------------------------------------

p_sweep <- ggplot(recurrence_sweep,
                  aes(x = min_tumor_samples, y = n_partners_retained)) +
  # Highlight bar at chosen cutoff
  geom_col(
    data = recurrence_sweep %>% filter(is_chosen),
    aes(x = min_tumor_samples, y = n_partners_retained),
    fill = "#c0392b", alpha = 0.15, width = 0.8
  ) +
  geom_vline(xintercept = recurrence_cutoff,
             linetype = "dashed", colour = "#c0392b", linewidth = 0.6) +
  annotate("text",
           x = recurrence_cutoff + 0.15, y = Inf,
           vjust = 1.6, hjust = 0,
           label = paste0("chosen cutoff: n >= ", recurrence_cutoff),
           size = 3.2, colour = "#c0392b") +
  geom_line(colour = "grey30", linewidth = 0.9) +
  geom_point(aes(colour = is_chosen), size = 3) +
  geom_text(
    aes(label = n_partners_retained),
    vjust = -0.8, size = 3.2, colour = "grey25"
  ) +
  scale_colour_manual(values = c("FALSE" = "grey30", "TRUE" = "#c0392b"),
                      guide = "none") +
  scale_x_continuous(
    breaks = sweep_thresholds,
    sec.axis = sec_axis(
      ~ . / n_tumor_total * 100,
      name = "% of tumor samples",
      breaks = sweep_thresholds / n_tumor_total * 100,
      labels = paste0(round(sweep_thresholds / n_tumor_total * 100), "%")
    )
  ) +
  scale_y_continuous(
    limits = c(0, max(recurrence_sweep$n_partners_retained) * 1.15),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title    = "Tumor-exclusive fusion partners retained at each recurrence cutoff",
    subtitle = paste0(
      "Tumor-exclusive pool: ",
      format(n_tumor_exclusive, big.mark = ","),
      " partners (n_normal = 0)  |  ",
      "Chosen cutoff n ≥ ", recurrence_cutoff,
      " retains ", format(n_filtered, big.mark = ","), " partners (",
      round(100 * n_filtered / n_tumor_exclusive, 1), "%)"
    ),
    x = "Minimum number of tumor samples partner must appear in (≥ n)",
    y = "Number of tumor-exclusive partners retained"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor    = element_blank(),
    plot.title.position = "plot"
  )

ggsave(
  file.path(arriba_dir, "tumor_exclusive_recurrence_sweep.png"),
  p_sweep,
  width = 10, height = 5.5, dpi = 300
)

message("Sweep plot written to: ", arriba_dir)


# ---- 11. Write outputs -------------------------------------------------------

write_csv(arriba_events,
          file.path(arriba_dir, "arriba_events_long.csv"))
write_csv(partner_long,
          file.path(arriba_dir, "fusion_partner_events_long.csv"))
write_csv(partner_summary,
          file.path(arriba_dir, "fusion_partner_summary_full.csv"))
write_csv(partner_tumor_exclusive,
          file.path(arriba_dir, "fusion_partner_tumor_exclusive.csv"))
write_csv(partner_summary_filtered,
          file.path(arriba_dir, "fusion_partner_summary_filtered.csv"))
write_csv(recurrence_sweep,
          file.path(arriba_dir, "tumor_exclusive_recurrence_sweep.csv"))

message("\n--- Section 3.4 complete ---")
message("Total unique partners (all tissues):         ", nrow(partner_summary))
message("Tumor-exclusive partners (n_normal = 0):     ", n_tumor_exclusive)
message("Partners after recurrence cutoff (n >= ", recurrence_cutoff, "): ",
        n_filtered,
        " (", round(100 * n_filtered / n_tumor_exclusive, 1), "%)")
message("Outputs written to: ", arriba_dir)