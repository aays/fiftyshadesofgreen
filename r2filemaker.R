#!usr/bin/env Rscript

# generates csv files out of plink table-formatted r2 output to then be used as input for zerocorrector.py
# usage - Rscript r2filemaker.R [filename]
# reference - notebook 8.7

library(tidyr)
library(dplyr)

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

outname <- paste(infile, '2', sep = '')

write.csv(chrdf, file = outname) 
