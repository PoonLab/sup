library(dplyr)
library(ggplot2)
library(ggrepel)

# read in sampled lineages
pang <- read.csv("data/pangolineages/ERR4363387_pangolineages.csv")
pang[1, ] # The actual sequence (not sampled)
# Compare to population
alldata <- read.csv("data/unc_covid/metadata.tsv", sep = "\t")

# Show the amount of diversity in lineage assignments
barplot(sort(table(pang$lineage), decreasing = TRUE))

# Show the most common lineages
sort(table(pang$lineage), decreasing = TRUE)[1:10]

# Compare to population ----
pangc <- table(pang$lineage)
pangdf <- data.frame(count = as.vector(pangc) / sum(pangc),
    lineage = names(pangc), df = "pang")
pangdf <- pangdf[pangdf$count > 0.005, ]

allc <- table(alldata$pangolin_lineage)
alldf <- data.frame(count = as.vector(allc) / sum(allc),
    lineage = names(allc), df = "all")
alldf <- alldf[alldf$lineage %in% pangdf$lineage, ]

both <- rbind(pangdf, alldf)

ggplot(both, aes(x = lineage, y = count, colour = df)) +
    geom_point() +
    labs(title = "Samples versus all lineages")
# Conclusion: NOT guessing randomly


# Better vis: incorporate boot prob somehow
head(pang)
pang2 <- lapply(unique(pang$lineage), function(x) {
    data.frame(lineage = x, prop = mean(pang$lineage == x))
}) %>%
    bind_rows() %>%
    right_join(pang, by = "lineage")
#pang2
ggplot(pang2) +
    aes(x = prop, y = probability,
        colour = lineage, label = lineage) +
    geom_point() +
    #geom_text_repel() +
    theme(legend.position = "none")

pangtab <- pang2 %>%
    group_by(prop, lineage) %>%
    summarise(y = 1, count = n(), .groups = "drop") %>%
    filter(prop > 0.025)

pang2 %>%
    # round to nearest 0.5
    #mutate(prop = round(prop*2, 1)/2,
    #    probability = round(probability*2, 1)/2) %>%
    group_by(prop, probability, lineage) %>%
    summarise(count = n(), .group = "drop") %>%
    ggplot() + theme_bw() +
        aes(x = prop, y = probability, colour = lineage, label = count) +
        geom_text() +
        theme(legend.position = "none") +
        annotate("text", x = pangtab$prop, y = 1, label = pangtab$lineage,
            hjust = 0.5, vjust = -1) +
        labs(x = "Proportion of Lineage in 10,00 Resampled Pangolin Outputs",
            y = "Pangolin's Reported Bootstrap Probability",
            title = "Pangolin Probabilities versus Resampled Lineage Calls",
            subtitle = "Numbers represent number of overlapping points;
                colours represent different lineages
                \nEach lineage will have a single proportion") +
        scale_x_continuous(breaks = seq(0, 1, 0.05))
