#!/usr/bin/env Rscript

# plots output from LDHelmet's rjmcmc + post_to_text steps

# usage: Rscript chromplot.R [input file (txt)] [outfile name] [outfile wd]
# e.g. Rscript chromplot.R chromosome_1.txt chromosome_1 analysis/graphs

# reference - notebook 5.3a

args <- commandArgs(trailingOnly = TRUE)

library(ggplot2)

# check for args
if (length(args) != 3) {
  stop("Missing an arg", call. = FALSE)
} 

# save args locally
infile <- read.csv(args[1], sep = ' ', header = TRUE, skip = 2)
outfile <- args[2]
graphwd <- args[3]

# fixes some format quirks with ldhelmet txt output
colfixer <- function(df) {
    colnames(df) <- colnames(df)[2:6]
    df <- df[,-6]
    return(df)
}

infile <- colfixer(infile)

# assert that infile is still a df
if (is.data.frame(infile) == FALSE) {
    stop("Error in reading inputted file", call. = TRUE) 
}

# create plot
plot <- ggplot(infile, aes(x = left_snp, y = mean)) +
  ggtitle(args[1]) +
  geom_line() +
  geom_smooth() +
  theme_bw()

outfile <- paste(outfile, '.png', sep = '')

ggsave(filename = outfile, plot = plot,
       path = graphwd,
    width = par("din")[1] * 1.75, height = par("din")[2],
    units = "in", dpi = 1000)   

