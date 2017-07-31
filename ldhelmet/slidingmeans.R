# does what slidingmeans.py does, but in R.
# takes in LDhelmet output and a windowsize
# ie. Rscript slidingmeans.R [ldhelmet file] [windowsize] [chromosome name]
# saves to 'genomedist.txt' in same directory (for now)
# for multiple chromosomes - for i in {1..n}; do Rscript slidingmeans.R chromosome_$i.txt 200000 chromosome_$i; done

# AH - 07/2017

library(dplyr, warn.conflicts = FALSE)
library(magrittr, warn.conflicts = FALSE)

args <- commandArgs(trailingOnly = TRUE)

filename <- as.character(args[1])
window <- as.integer(args[2])
chrname <- as.character(args[3])

windowcalc <- function(df, windowsize) {
    df %<>% select(left_snp = V1, right_snp = V2, rho = V3) %>% # get cols
        mutate(diff = right_snp - left_snp) %>% # get diff
        mutate(scaledrho = rho * diff) %>% # use diff for scaled rho
        mutate(window = ceiling(left_snp/windowsize)) %>% # compute windows
        group_by(window) %>% # group by windows
        summarise(rho = mean(scaledrho)) %>% # get means
        mutate(chr = chrname) %>%
        select(-window)
    return(df)
}

df <- read.csv(filename, sep = ' ', skip = 3, header = FALSE)
df %<>% windowcalc(window)

df <- df[,c(2,1)] # reorder

names(df) <- NULL # suppress colnames in printing
write.table(df, file = 'genomedist.txt', sep = ' ', 
            col.names = FALSE, quote = FALSE, append = TRUE)
