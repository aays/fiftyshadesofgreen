#!/usr/bin/env Rscript

# usage - Rscript r2graphweirhillv2.R [input file] [outfile name] [plot wd]

# reference - notebook 8.7

library(ggplot2)
library(tidyr)
library(dplyr)

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
    weirhilleq <- '((10 + p*d)/(22 + (13*p*d) + (p*d)^2))*(1 + (((3 + (p*d))/(24*(22 + (13*p*d) + (p*d)^2)^2)) * (12 + (12*p*d) + (p*d)^2)))'

    subdf <- df %>%
        select(one_of(c('d', 'r2'))) %>%
        apply(2, as.numeric)
    subdf <- as.data.frame(subdf)

    tempdf <- as.data.frame(predict(nls(paste('r2', '~', weirhilleq),
                      data = subdf, control = list(maxiter = 500), start = list(p = 0.5))))
    tempdf$d <- df$d
    colnames(tempdf)[1] <- 'rho'
    tempdf <- apply(tempdf, 2, as.numeric)

    r2plot <- ggplot(subdf, aes(x = d, y = r2)) +
        geom_point() +
        geom_line(data = as.data.frame(tempdf), aes(x = d, y = rho), col = 'red') +
        geom_smooth()

    ggsave(filename = paste(outfile,'zeroed.png', sep = ''),
           plot = r2plot, path = graphwd,
           width = par("din")[2], height = par("din")[2],
           units = "in", dpi = 750)

    paste(outfile, summary(nls(paste('r2', '~', weirhilleq),
                      data = subdf, control = list(maxiter = 500),
                      start = list(p = 0.5)))$p)
}

data <- read.csv(infile, header = TRUE)

colnames(data) <- 'column'

data <- mutate(data, column = gsub(pattern = "[ ]{1,}", replacement = ",", x = column)) %>%
separate(col = "column", into = c('pos1', 'pos2', 'r2'), sep = ",")

chrdf <- diffgetter(data)

plotter(chrdf, outfile, graphwd)

