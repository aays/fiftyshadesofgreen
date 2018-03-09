#!usr/bin/env Rscript

# generates csv files out of plink table-formatted r2 output to then be used as input for zerocorrector.py

# usage - Rscript r2corrector.R [filename]
# eg r2corrector.R chromsome1.ld will generate a file called chromosome1.ldz in the same dir

# AH - 04/2017

library(tidyr)
library(dplyr, warn.conflicts = FALSE)

args <- commandArgs(trailingOnly = TRUE)

# check for args
if (length(args) != 1) {
  stop("Missing an arg", call. = FALSE)
} 

infile <- args[1]

remover <- function(df){
  df <- df[,c(2,5,7)]
  return(df)
} 

data <- read.csv(infile, header = TRUE) 
colnames(data) <- 'column'
data <- mutate(data, column = gsub(pattern = "[ ]{1,}", replacement = ",", x = column)) %>%
separate(col = "column", into = c('chr1', 'pos1', '.', 'chr2', 'pos2', '.', 'r2'), sep = ",") 

chrdf <- remover(data)

outname <- paste(infile, 'z', sep = '')

write.csv(chrdf, file = outname) 
