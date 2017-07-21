#!/usr/bin/env Rscript

# functions for taking in table-formatted output from plink's --r2 function and generates a line
# that fit an equation modelling LD decay over distance (Weir & Hill 1986, Cutter et al. 2006)

# assumes n = 24

# AH - 07/2017

library(tidyr)
library(dplyr, warn.conflicts = FALSE)

diffgetter <- function(df) {
    colnames(df) <- 'column'
    df <- df %>% 
      mutate(column = gsub(pattern = "[ ]{1,}", replacement = ",", x = column)) %>%
      separate(col = "column", into = c('pos1', 'pos2', 'r2'), sep = ",")
      mutate(d = as.numeric(pos2) - as.numeric(pos1)) %>%
      select(one_of(c('d', 'r2'))) %>%
      apply(2, as.numeric)
    
    df[,c(1,2)] <- sapply(df[,c(1,2)], as.numeric)
    return(as.data.frame(df))
}

linemaker <- function(df) {
    df <- as.data.frame(df)
    weirhilleq <- '((10 + p*d)/(22 + (13*p*d) + (p*d)^2))*(1 + (((3 + (p*d))/(24*(22 + (13*p*d) + (p*d)^2)^2)) * (12 + (12*p*d) + (p*d)^2)))'

    subdf <- df %>%
        select(one_of(c('d', 'r2'))) %>%
        apply(2, as.numeric) %>%
        as.data.frame()

    linedf <- nls(paste('r2', '~', weirhilleq),
                  data = subdf, control = list(maxiter = 500), start = list(p = 0.5)) %>%
              predict() %>%
              as.data.frame()
  
    message('Model fit complete.')
    
    linedf$d <- df$d
    colnames(linedf)[1] <- 'rho'
    linedf <- apply(tempdf, 2, as.numeric)
    
    paste(outfile, summary(nls(paste('r2', '~', weirhilleq),
                      data = subdf, control = list(maxiter = 500),
                      start = list(p = 0.5)))$p)
    
    return(as.data.frame(linedf))
                      
}
