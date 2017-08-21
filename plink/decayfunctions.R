# useful functions for working with LD decay data

# AH - 08/2017

decaypoint <- function(df, frac, verbose = FALSE) {
    weirhilleq <- '((10 + p*d)/(22 + (13*p*d) + (p*d)^2))*(1 + (((3 + (p*d))/(24*(22 + (13*p*d) + (p*d)^2))) * (12 + (12*p*d) + (p*d)^2)))'

    subdf <- df %>%
        mutate(d = pos2 - pos1) %>%
        select(d, r2) %>%
        apply(2, as.numeric) %>%
        as.data.frame()

    linedf <- nls(paste('r2', '~', weirhilleq),
                  data = subdf, control = list(maxiter = 500), start = list(p = 0.5)) %>%
              predict() %>%
              as.data.frame()
  
    linedf$d <- subdf$d
    colnames(linedf)[1] <- 'rho'
    
    start.rho <- linedf$rho[1]
    start.d <- linedf$d[1]
    stop.rho <- start.rho * frac
    
    linedf %<>% arrange(desc(rho))
    
    snippet <- filter(linedf, rho >= stop.rho)
    
    end.rho <- snippet$rho[dim(snippet)[1]]
    end.d <- snippet$d[dim(snippet)[1]]
    
    if (verbose == FALSE) {
        out <- c(start.rho, end.rho, start.d, end.d, end.d - start.d)
    } else if (verbose == TRUE) {
        out <- paste('start rho = ', start.rho,
                    '\nend.rho = ', end.rho,
                    '\nstart.d = ', start.d,
                    '\nend.d = ', end.d,
                    '\ndiff =', end.d - start.d)
    }
    
    return(out)
}

getline <- function(df) {
    # provides decay line as a df
    weirhilleq <- '((10 + p*d)/(22 + (13*p*d) + (p*d)^2))*(1 + (((3 + (p*d))/(24*(22 + (13*p*d) + (p*d)^2))) * (12 + (12*p*d) + (p*d)^2)))'

    subdf <- df %>%
        mutate(d = pos2 - pos1) %>%
        select(d, r2) %>%
        apply(2, as.numeric) %>%
        as.data.frame()

    linedf <- nls(paste('r2', '~', weirhilleq),
                  data = subdf, control = list(maxiter = 500), start = list(p = 0.5)) %>%
              predict() %>%
              as.data.frame()
  
    linedf$d <- subdf$d
    colnames(linedf)[1] <- 'rho'
    linedf %<>% arrange(desc(rho))
    
    return(linedf)
}

decaypoint.line <- function(linedf, frac, verbose = FALSE) {
    # takes input from getline() instead of direct plink data
    start.rho <- linedf$rho[1]
    start.d <- linedf$d[1]
    stop.rho <- start.rho * frac
    
    linedf %<>% arrange(desc(rho))
    
    snippet <- filter(linedf, rho >= stop.rho)
    
    end.rho <- snippet$rho[dim(snippet)[1]]
    end.d <- snippet$d[dim(snippet)[1]]
    
    if (verbose == FALSE) {
        out <- c(start.rho, end.rho, start.d, end.d, end.d - start.d)
    } else if (verbose == TRUE) {
        out <- paste('start rho = ', start.rho,
                    '\nend.rho = ', end.rho,
                    '\nstart.d = ', start.d,
                    '\nend.d = ', end.d,
                    '\ndiff =', end.d - start.d)
    }
    
    return(out)
}

decaypoint.to.df <- function(linedf, frac, chrname) {
    # takes input from decayline and returns in df format
    start.rho <- linedf$rho[1]
    start.d <- linedf$d[1]
    stop.rho <- start.rho * frac
    
    linedf %<>% arrange(desc(rho))
    
    snippet <- filter(linedf, rho >= stop.rho)
    
    end.rho <- snippet$rho[dim(snippet)[1]]
    end.d <- snippet$d[dim(snippet)[1]]
    
    out <- data.frame(chr = chrname, 
                     rho.start = start.rho,
                     rho.end = end.rho,
                     d.start = start.d,
                     d.end = end.d,
                     diff = end.d - start.d)
    return(out)
}

# halfdecay <- getline(dflist[[1]]) %>% decaypoint.to.df(0.5, names(dflist)[[1]])
# for (i in c(2:17)){ 
#     tempdf <- getline(dflist[[i]]) %>% decaypoint.to.df(0.5, names(dflist)[[i]])
#     halfdecay %<>% bind_rows(tempdf)
# }

flatline.point <- function(linedf, chrname) {
    # uses rates of change to approximately determine when decay 'flatlines'
    linedf %<>%
        group_by(d) %>%
        summarise(mean.rho = mean(rho)) %>%
        mutate(d2 = lead(d), mean.rho2 = lead(mean.rho)) %>%
        mutate(ddx = abs(round((mean.rho2 - mean.rho)/(d2 - d), digits = 3))) # approx rate of change
    
    start.rho <- linedf$mean.rho[1]
    start.d <- linedf$d[1]
    
    snippet <- filter(linedf, ddx > 0)
    
    end.rho <- snippet$mean.rho[dim(snippet)[1]]
    end.d <- snippet$d[dim(snippet)[1]]
    
    out <- data.frame(chr = chrname, 
                     rho.start = start.rho,
                     rho.end = end.rho,
                     d.start = start.d,
                     d.end = end.d,
                     diff = end.d - start.d)
    return(out)
}
