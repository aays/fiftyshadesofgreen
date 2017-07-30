# Primes LDhelmet output for use with scaledmeans.py
# takes in an LDhelmet outfile and scales each rho value to the distance between SNPs

# AH - 07/2017

library(dplyr, warn.conflicts = FALSE)
library(magrittr, warn.conflicts = FALSE)

args <- commandArgs(trailingOnly = TRUE)

infile <- args[1]

df <- read.csv(infile, sep = ' ', skip = 3, header = FALSE)

df %<>% 
    select(left_snp = V1, right_snp = V2, mean = V3) %>% 
    mutate(diff = right_snp - left_snp) %>% 
    mutate(scaledmean = diff * mean) %>%
    select(-mean)

write.table(df, file = paste('scaled', infile, sep =''), sep = ' ')
