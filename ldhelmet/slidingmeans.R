# does what slidingmeans.py does, but in R.
# takes in LDhelmet output and a windowsize
# ie. Rscript slidingmeans.R [ldhelmet file] [windowsize]
# saves to same directory (for now)

# AH - 07/2017

library(dplyr, warn.conflicts = FALSE)
library(magrittr, warn.conflicts = FALSE)

args <- commandArgs(trailingOnly = TRUE)

filename <- as.character(args[1])
window <- as.integer(args[2])

windowcalc <- function(df, windowsize) {
    df %<>% select(left_snp = V1, right_snp = V2, rho = V3) %>% # get cols
        mutate(diff = right_snp - left_snp) %>% # get diff
        mutate(scaledrho = rho * diff) %>% # use diff for scaled rho
        mutate(window = ceiling(left_snp/windowsize)) %>% # compute windows
        group_by(window) %>% # group by windows
        summarise(rho = mean(scaledrho)) # get means
    return(df)
}

df <- read.csv(filename, sep = ' ', skip = 3, header = FALSE)
df %<>% windowcalc(window)

write.table(df, file = paste('means_', filename, sep = ''), sep = ' ')  
