library(dplyr)
library(tidyr)
library(ggplot2)
theme_set(theme_bw())

# Cutoff - anything above this is basically certain
cutoff <- 0.99

# Load data
#thissam <- readRDS("data/unc_covid/open_seSAMe-S-ERR4364007-1.RDS")
uncsam <- list.files("data/unc_covid", pattern = "*\\.RDS")

plotlist <- vector(mode = "list", length = length(uncsam))
meaningful <- c()
for(i in 1:length(uncsam)){
    thissam <- readRDS(paste0("data/unc_covid/", uncsam[[i]]))
    if(is.null(dim(thissam))) next

    # My algorithm doesn't actually assign probabilities to gaps
    thissam <- thissam[,1:4]

    # Find the probability of the consensus
    maxprob <- apply(thissam, 1, max)
    sumprob <- apply(thissam, 1, sum)
    thissam <- apply(thissam, 2, function(x) x / (sumprob + 0.00001))


    samdf <- data.frame(thissam) %>% 
        mutate(id = 1:n()) %>%
        slice(which(maxprob < cutoff)) %>%
        pivot_longer(cols = -id, names_to = "alph", values_to = "q") %>%
        filter(q > 0)
    
    meaningful[i] <- length(unique(samdf$q)) > 2
    plotlist[[i]] <- ggplot(samdf) + 
        aes(x = factor(id), y = alph, fill = q) + 
        geom_tile() + 
        scale_fill_viridis_c(option = "A", direction = -1) +
        labs(title = uncsam[[i]], x = "Locus", y = NULL)
}

cowplot::plot_grid(plotlist = plotlist[meaningful], ncol = 2)
