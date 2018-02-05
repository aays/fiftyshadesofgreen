#!/usr/bin/env Rscript

# takes in table-formatted output from plink's --r2 function and generates plots
# that fit an equation modelling LD decay over distance (Weir & Hill 1986, Cutter et al. 2006)

# assumes n = 24

# usage - Rscript lddecaywh.R [input file] [outfile name] [plot wd]
# eg Rscript lddecaywh.R chrom1.txt chrom1 /analysis/graphs creates chrom1.png in /analysis/graphs

# AH - 03/2017

# update 02/2018 - fixed equation

library(ggplot2)
library(tidyr)
library(dplyr, warn.conflicts = FALSE)

args <- commandArgs(trailingOnly = TRUE)

# check for args
if (length(args) != 3) {
  stop("Missing an arg", call. = FALSE)
}

infile <- args[1]
outfile <- args[2]
graphwd <- args[3]

diffgetter <- function(df) {
    df <- df %>%
      mutate(d = as.numeric(pos2) - as.numeric(pos1)) %>%
      select(one_of(c('d', 'r2'))) %>%
      apply(2, as.numeric)
    df[,c(1,2)] <- sapply(df[,c(1,2)], as.numeric)
    return(as.data.frame(df))
}

plotter <- function(df, outfile, graphwd) {
    df <- as.data.frame(df)
    weirhilleq <- '((10 + p*d)/(22 + (13*p*d) + (p*d)^2))*(1 + (((3 + (p*d))/(24*(22 + (13*p*d) + (p*d)^2))) * (12 + (12*p*d) + (p*d)^2)))'

    subdf <- df %>%
        select(one_of(c('d', 'r2'))) %>%
        apply(2, as.numeric) %>%
        as.data.frame()

    tempdf <- nls(paste('r2', '~', weirhilleq),
                  data = subdf, control = list(maxiter = 500), start = list(p = 0.5)) %>%
              predict() %>%
              as.data.frame()
  
    message('Model fit complete.')
    
    tempdf$d <- df$d
    colnames(tempdf)[1] <- 'rho'
    tempdf <- apply(tempdf, 2, as.numeric)
    
    #message('Creating plot...')
    
    #r2plot <- ggplot(subdf, aes(x = d, y = r2)) +
    #    geom_point() +
    #    geom_line(data = as.data.frame(tempdf), aes(x = d, y = rho), col = 'red')

    #message(paste('Saving plot to ', graphwd, outfile, '...', sep = ''))
  
    #ggsave(filename = paste(outfile, '.png', sep = ''),
    #       plot = r2plot, path = graphwd,
    #       width = par("din")[2], height = par("din")[2],
    #       units = "in", dpi = 750)
    
    message('Complete.')
    
    paste(outfile, summary(nls(paste('r2', '~', weirhilleq),
                      data = subdf, control = list(maxiter = 500),
                      start = list(p = 0.5)))$p)
}

message('Loading in file...')

data <- read.csv(infile, header = TRUE)

colnames(data) <- 'column'

message('Complete.')

data <- mutate(data, column = gsub(pattern = "[ ]{1,}", replacement = ",", x = column)) %>%
separate(col = "column", into = c('pos1', 'pos2', 'r2'), sep = ",")

chrdf <- diffgetter(data)

plotter(chrdf, outfile, graphwd)

