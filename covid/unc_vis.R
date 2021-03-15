library(dplyr)
library(here)
library(tidyr)
library(ggplot2)
library(ggthemes)
theme_set(theme_bw())

# Cutoff - anything above this is basically certain
cutoff <- 0.999

# Load data
uncsam <- list.files(here("data", "unc_covid"), pattern = "*\\.RDS")

plotlist <- vector(mode = "list", length = length(uncsam))
meaningful <- c()
for (i in seq_along(uncsam)) {
    thissam <- readRDS(here("data", "unc_covid", uncsam[[i]]))
    if (is.null(dim(thissam))) next

    # My algorithm doesn't actually assign probabilities to gaps
    thissam <- thissam[, 1:4]

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
        scale_fill_gradient2(low = "#276419", mid = "#f7f7f7", 
            high = "#8e0152", midpoint = 0.5,
            limits = c(0, 1)) +
        labs(title = uncsam[[i]], x = NULL, y = NULL) +
        theme(axis.text.x = element_text(angle = 90, size = 7))
}

meaningful[is.na(meaningful)] <- FALSE
print(patchwork::wrap_plots(plotlist[meaningful], 
    ncol = 2, guides = "collect"))
