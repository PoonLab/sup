library(tidyr)
library(dplyr)
library(ggplot2) ; theme_set(theme_bw())

source('utils-zanini.R')


# ---- Data ----

pat <- 1
tpt <- c(1:5,0)
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


# ---- Estimation ----


# Select only the positions that have
# a probability above the threshold:
thresh <- 0.999
dft <- df %>%
    mutate(mono.a = A>thresh,
           mono.c = C>thresh,
           mono.g = G>thresh,
           mono.t = T>thresh)


# variance of presence across time points:
dfv <- dft %>%
    group_by(pos) %>%
    summarise(cv.a = sd(A)/mean(A),
              cv.c = sd(C)/mean(C),
              cv.g = sd(G)/mean(G),
              cv.t = sd(T)/mean(T),
              n=n()) %>%
    mutate()
dfv
    
mean(dfs$is.invar)
