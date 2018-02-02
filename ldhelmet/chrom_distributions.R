# takes in files outputted from find_hotspots.py (Singhal 2015)

library(dplyr, warn.conflicts = FALSE)
library(magrittr)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)

s_file <- args[1]

d <- read.csv(s_file)

chrnames <- unique(d$chr)

for (chrom in chrnames){

    plot <- d %>%
        filter(chr == chrom) %>%
        select(block_rate) %>%
        ggplot(aes(x = block_rate, y = ..count../sum(..count..))) +
        geom_histogram(binwidth = 0.01) +
        labs(x = 'recombination rate', y = 'proportion') +
        theme_bw() +
        geom_vline(aes(xintercept = 5 * 0.00443298965947912), linetype = 'dashed', col = 'red') +
        geom_vline(aes(xintercept = 0.00443298965947912), linetype = 'dashed', col = 'blue') +
        ggtitle(paste(chrom, 'binwidth = 0.01'))
    
    ggsave(paste(chrom, '.png'), plot = plot, path = getwd(), dpi = 300)
    message(paste('done', chrom))
}
