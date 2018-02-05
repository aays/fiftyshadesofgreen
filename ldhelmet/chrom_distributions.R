# takes in files outputted from find_hotspots.py (Singhal 2015)

library(dplyr, warn.conflicts = FALSE)
library(magrittr)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)

s_file <- args[1]

d <- read.csv(s_file)

chrnames <- unique(d$chr)

refdf <- data.frame(chromosome = c('chromosome_1', 'chromosome_10', 'chromosome_11', 'chromosome_12', 
                                   'chromosome_13', 'chromosome_14', 'chromosome_15', 'chromosome_16', 
                                   'chromosome_17', 'chromosome_2', 'chromosome_3', 'chromosome_4', 'chromosome_5', 
                                   'chromosome_6', 'chromosome_7', 'chromosome_8', 'chromosome_9'),
                    rho = c(0.004463396, 0.004497173, 0.005309452, 0.002591899, 0.005389904, 0.008493510, 
                            0.011547150, 0.004134429, 0.006021989, 0.004205883, 0.003980304, 0.003879390, 
                            0.005102792, 0.003220011, 0.002944370, 0.004774266, 0.003762676),
                    lengths = c(8033585, 6576019, 3826814, 9730733, 5206065, 4157777, 1922860, 7783580, 
                                7188315, 9223677, 9219486, 4091191, 3500558, 9023763, 6421821, 5033832, 7956127))

rhogetter <- function(chr, indf){
    out <- indf %>%
        filter(chromosome == chr) %>%
        select(rho) %>%
        unlist()
    return(out)
}


for (chrom in chrnames){

    plot <- d %>%
        filter(chr == chrom) %>%
        select(block_rate) %>%
        ggplot(aes(x = block_rate, y = ..count../sum(..count..))) +
        geom_histogram(binwidth = 0.01) +
        labs(x = 'recombination rate', y = 'proportion') +
        theme_bw() +
        geom_vline(aes(xintercept = 5 * rhogetter(chrom, refdf)), linetype = 'dashed', col = 'red') +
        geom_vline(aes(xintercept = rhogetter(chrom, refdf)), linetype = 'dashed', col = 'blue') +
        ggtitle(paste(chrom, 'binwidth = 0.01'))
    
    ggsave(paste(chrom, '.png'), plot = plot, path = getwd(), dpi = 300)
    message(paste('done', chrom))
}
