
```
library(regioneR)
library(GenomicRanges)
library(rtracklayer)
library(ggplot2)
library(dplyr)
library(ChIPseeker)
library(txdbmaker)
library(GenomicFeatures)
library(tidyverse)

setwd("~/dissertation/RIANseq_files")

peak_files <- list.files("narrowPeaks", pattern = "_peaks.narrowPeak$", full.names = TRUE)
sample_names <- gsub("_peaks.narrowPeak", "", basename(peak_files))

peak_list <- setNames(
  lapply(peak_files, function(f) toGRanges(f)),
  sample_names
)

dsb <- toGRanges("DSB_coordinates_CHM13.bed")
fai <- read.table("chm13v2.0.fa.fai", header = FALSE)
genome_df <- data.frame(chr = fai$V1, start = 1, end = fai$V2)
genome <- toGRanges(genome_df)

sapply(peak_list, length)

lapply(peak_list, function(p) unique(as.character(seqnames(p))))

sapply(peak_list, function(p) range(width(p)))

sapply(peak_list, function(p) any(is.na(start(p)) | is.na(end(p))))

chrM_summary <- data.frame(
  sample = sample_names,
  total_peaks = sapply(peak_list, length),
  chrM_peaks = sapply(peak_list, function(p) sum(as.character(seqnames(p)) == "chrM"))
) %>% mutate(chrM_pct = round(100 * chrM_peaks / total_peaks, 1))

print(chrM_summary)

peak_list_clean <- lapply(peak_list, function(p) p[seqnames(p) != "chrM"])

hsat2_region <- GRanges("chr1", IRanges(127419673, 127420254))
sapply(peak_list_clean, function(p) sum(countOverlaps(p, hsat2_region) > 0))

flagstat_files <- list.files("flagstat_rianseq", pattern = "\\.flagstat\\.txt$", full.names = TRUE)

parse_flagstat <- function(path) {
  lines <- readLines(path)
  
  mapped_line <- grep("mapped \\(", lines, value = TRUE)[1]
  mapped_pct <- as.numeric(sub(".*\\(([0-9.]+)%.*", "\\1", mapped_line))
  
  paired_line <- grep("properly paired", lines, value = TRUE)[1]
  paired_pct <- as.numeric(sub(".*\\(([0-9.]+)%.*", "\\1", paired_line))
  
  total_reads <- as.numeric(sub("^([0-9]+).*", "\\1", lines[1]))
  
  data.frame(
    mapped_pct = mapped_pct,
    properly_paired_pct = paired_pct,
    total_reads = total_reads
  )
}

flagstat_table <- do.call(rbind, lapply(flagstat_files, function(f) {
  sample <- gsub("\\.flagstat\\.txt$", "", basename(f))
  cbind(sample = sample, parse_flagstat(f))
}))

print(flagstat_table)

insert_sizes <- read.delim("rianseq_stats/insert_size_summary.txt") %>%
  rename(sample = Sample)

peak_counts <- read.delim(
  "rianseq_results/summary/peak_counts_summary.txt",
  header = FALSE,
  col.names = c("sample", "total_peaks")
)

peak_counts <- peak_counts %>%
  mutate(sample = gsub("_peaks$", "", sample))

table_3_4 <- insert_sizes %>%
  left_join(peak_counts, by = "sample") %>%
  left_join(chrM_summary, by = "sample") %>%
  left_join(flagstat_table, by = "sample")

chrM_summary_join <- chrM_summary %>%
  dplyr::select(sample, chrM_peaks, chrM_pct)

table_3_4 <- insert_sizes %>%
  left_join(peak_counts, by = "sample") %>%
  left_join(chrM_summary_join, by = "sample") %>%
  left_join(flagstat_table, by = "sample")

write.csv(table_3_4, "table_3_4_alignment_qc.csv", row.names = FALSE)

flagged <- c(
  "diva_siCtrl_DMSO",
  "diva_siUpf1_tam",
  "rpe1_upf1ko_etoposide"
)

table_3_5 <- table_3_4 %>%
  mutate(group = ifelse(sample %in% flagged, "flagged", "other")) %>%
  group_by(group) %>%
  summarise(
    n = n(),
    mean_peaks = mean(total_peaks),
    mean_chrM_pct = mean(chrM_pct)
  )

print(table_3_5)

get_fold_enrichment_range <- function(xls_path) {
  df <- read.delim(xls_path, comment.char = "#")
  range(df$fold_enrichment, na.rm = TRUE)
}

flagged <- c(
  "diva_siCtrl_DMSO",
  "diva_siUpf1_tam",
  "rpe1_upf1ko_etoposide"
)

fold_enrichment_ranges <- sapply(flagged, function(s) {
  path <- file.path("rianseq_peaks", paste0(s, "_peaks.xls"))
  rng <- get_fold_enrichment_range(path)
  paste0(round(rng[1], 1), "-", round(rng[2], 1))
})

print(fold_enrichment_ranges)

fig_data <- table_3_4 %>%
  mutate(flagged = sample %in% flagged)

ggplot(fig_data, aes(x = reorder(sample, total_peaks), y = total_peaks, fill = flagged)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("FALSE" = "grey60", "TRUE" = "firebrick")) +
  labs(x = NULL, y = "Total peaks (chrM excluded)", fill = "Flagged sample") +
  theme_minimal()

ggsave("figures/figure_3_5_peak_counts.png", width = 7, height = 5)

head(peak_counts)

table_3_6 <- data.frame(
  sample = c(
    "diva_siCtrl_DMSO",
    "diva_siCtrl_DMSO_RNaseH",
    "diva_siCtrl_tam",
    "diva_siCtrl_tam_RNaseH",
    "diva_siUpf1_DMSO",
    "diva_siUpf1_tam",
    "diva_siUpf1_tam_RNaseH",
    "rpe1_WT_DMSO",
    "rpe1_WT_etoposide",
    "rpe1_WT_etoposide_RNaseH",
    "rpe1_upf1ko_DMSO",
    "rpe1_upf1ko_etoposide",
    "rpe1_upf1ko_etoposide_RNaseH"
  ),
  reads_in_peaks = c(
    18229, 29997, 95095, 29036, 56924, 15083, 164910,
    34877, 88602, 61606, 60562, 18744, 15189
  ),
  total_reads = c(
    6175238, 2645114, 1428948, 6815436, 3400932, 7597958, 5457902,
    6460436, 5304262, 5694588, 5860914, 6783266, 3521216
  ),
  pct = c(
    0.30, 1.13, 6.65, 0.43, 1.67, 0.20, 3.02,
    0.54, 1.67, 1.08, 1.03, 0.28, 0.43
  )
)

write.csv(table_3_6, "table_3_6_enrichment.csv", row.names = FALSE)

table_3_7 <- read.delim(
  "summary/yeast_spikein_summary.txt",
  header = FALSE,
  col.names = c(
    "sample",
    "yeast_raw",
    "yeast_filtered",
    "yeast_filtered_per_million"
  )
)

print(table_3_7)

diva_samples <- c(
  "diva_siCtrl_DMSO",
  "diva_siCtrl_tam",
  "diva_siCtrl_DMSO_RNaseH",
  "diva_siCtrl_tam_RNaseH",
  "diva_siUpf1_DMSO",
  "diva_siUpf1_tam",
  "diva_siUpf1_tam_RNaseH"
)

all_peak_chroms <- unique(unlist(lapply(
  peak_list_clean[diva_samples],
  function(p) as.character(seqnames(p))
)))

dsb_chroms <- unique(as.character(seqnames(dsb)))
shared_chroms <- union(all_peak_chroms, dsb_chroms)
genome_subset <- genome[genome_df$chr %in% shared_chroms, ]

results <- list()
failed_samples <- c()

for (s in diva_samples) {
  cat("Running:", s, "\n")
  
  pt <- tryCatch({
    permTest(
      A = dsb,
      B = peak_list_clean[[s]],
      ntimes = 1000,
      randomize.function = randomizeRegions,
      evaluate.function = meanDistance,
      genome = genome_subset
    )
  }, error = function(e) {
    cat("FAILED:", s, "-", conditionMessage(e), "\n")
    failed_samples <<- c(failed_samples, s)
    return(NULL)
  })
  
  if (!is.null(pt)) results[[s]] <- pt
}

saveRDS(results, "meanDistance_results.rds")

summary_table <- data.frame(
  sample = names(results),
  observed_mean_dist = sapply(results, function(x) x$meanDistance$observed),
  permuted_mean_dist = sapply(results, function(x) mean(x$meanDistance$permuted)),
  p_value = sapply(results, function(x) x$meanDistance$pval),
  z_score = sapply(results, function(x) x$meanDistance$zscore)
)

summary_table$p_adjusted <- p.adjust(summary_table$p_value, method = "BH")

write.csv(summary_table, "table_3_8_permutation.csv", row.names = FALSE)
print(summary_table)

genome_df$chr %in% shared_chroms

exact_overlap <- read.delim("rianseq_results/summary/diva_dsb_overlap_summary.txt")
proximity <- read.delim("rianseq_results/summary/diva_dsb_proximity_summary.txt")

table_3_8 <- exact_overlap %>%
  left_join(proximity, by = "Sample") %>%
  left_join(summary_table, by = c("Sample" = "sample"))

write.csv(table_3_8, "table_3_8_full.csv", row.names = FALSE)

diva_fig_data <- table_3_4 %>%
  filter(sample %in% diva_samples)

ggplot(diva_fig_data, aes(x = sample, y = total_peaks)) +
  geom_col() +
  coord_flip() +
  labs(x = NULL, y = "Total peaks (chrM excluded)") +
  theme_minimal()

ggsave("figure_3_6_diva_peaks.png", width = 7, height = 5)

plot(results[["diva_siCtrl_tam"]])

rpe1_summary <- read.delim("rianseq_results/summary/rpe1_genomewide_summary.txt")

write.csv(rpe1_summary, "table_3_9_rpe1.csv", row.names = FALSE)

ggplot(rpe1_summary, aes(x = Sample, y = Total_Peaks)) +
  geom_col() +
  coord_flip() +
  labs(x = NULL, y = "Total peaks") +
  theme_minimal()

ggsave("figure_3_8_rpe1_peaks.png", width = 7, height = 5)

txdb <- txdbmaker::makeTxDbFromGFF("chm13v2.0.gtf", format = "gtf")

anno_list <- lapply(peak_list_clean, function(p) {
  annotatePeak(
    p,
    TxDb = txdb,
    tssRegion = c(-3000, 3000),
    verbose = FALSE
  )
})

plotAnnoBar(anno_list)

ggsave("figure_annotation_rianseq.png", width = 9, height = 6)

anno_table <- lapply(names(anno_list), function(s) {
  df <- as.data.frame(anno_list[[s]]@annoStat)
  df$sample <- s
  df
}) %>% bind_rows()

write.csv(anno_table, "table_annotation_rianseq.csv", row.names = FALSE)
```