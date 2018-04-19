# Rscript make_hotspot_dist2k.R dist2k_full.txt
# quick and dirty script to make input for antr_diversity_hotspots.py

library(dplyr)
library(magrittr)
library(readr)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
file <- args[1]

d <- read_csv(file)

non_overlapping <- function(df, windowsize) {
    df %<>%
        mutate(div = floor(block_start / windowsize), div2 = lead(div)) %>%
        filter(div == div2) %>%
        select(-contains('div'))
    return(df)
}

d %<>%
    non_overlapping(2000) %>%
    select(-flank_rate, -rate_ratio) %>%
    group_by(chr) %>%
    mutate(chrom_mean = mean(block_rate, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(actual_ratio = block_rate / chrom_mean) %>%
    mutate(is_hotspot = case_when(
        actual_ratio >= 5.0 ~ 1,
        actual_ratio < 5.0 ~ 0))


outname <- str_replace(file, '.txt', '_hotspots.txt')
write.csv(d, outname, quote = FALSE, row.names = FALSE)
