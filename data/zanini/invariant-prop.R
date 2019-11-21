library(tidyr)
library(dplyr)
library(ggplot2) ; theme_set(theme_bw())

source('utils-zanini.R')


# ---- Data

#' Estimates the proportion of invariant position bases
#' for one patient from her/his longitudinal sequences (Zanini's data)
#' @param pat Integer. Patient number
#' @param tpt Integer vector. Time points of the sampled sequences to include.
#' @param thresh Numeric. Threshold frequency of a nucleotide to consider a position invariant.
estim_prop_invariant <- function(pat, tpt, thresh) { #pat=2
    
    print(paste('Patient',pat,'...'))
    
    nt <- length(tpt)
    dfl <- list()
    
    for (i in seq_along(tpt)) {
        dfl[[i]] <- read_data(patient = pat, 
                              timepoint = tpt[i])
        dfl[[i]]$tpt <- tpt[i]
        nr <- nrow(dfl[[i]])
        dfl[[i]]$pos <- 1:nr
    }
    
    # Convert to probabilities:
    df <- do.call('rbind', dfl) %>%
        select(-d, -N) %>%
        mutate(n = A+C+G+T)
    df[,1:4] <- df[,1:4] / df$n
    
    df <- select(df, -n)
    
    # ---- Estimation
    
    df2 <- df %>%
        pivot_longer(-c('tpt','pos'),
                     names_to = 'base',
                     values_to = 'proba') %>%
        mutate(test1 = proba > thresh)
    
    dfs <- df2 %>% 
        group_by(pos, base) %>%
        summarize(a = sum(test1)) %>%
        mutate(base.is.invar = (a==nt))
    
    df.inv <- dfs %>%
        group_by(pos) %>%
        summarize(pos.invar = max(base.is.invar))
    
    res <- mean(df.inv$pos.invar, na.rm = TRUE)
    return(res)
}

pat.vec <- c(1:11)
tpt <- c(1:2,0)
thresh <- 0.99


x <- sapply(pat.vec, estim_prop_invariant, tpt=tpt,thresh=thresh)

hist(x)
