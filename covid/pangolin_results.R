library(ggplot2)

# read in sampled lineages
pang <- read.csv("data/pangolineages/SRR13020989_pangolineages.csv")
pang[1,] # The actual sequence (not sampled)
# Compare to population
alldata <- read.csv("data/unc_covid/metadata.tsv", sep = "\t")

# Show the amount of diversity in lineage assignments
barplot(sort(table(pang$lineage), decreasing = TRUE))

# Show the most common lineages
sort(table(pang$lineage), decreasing = TRUE)[1:10]

# Compare to population ----
pangc <- table(pang$lineage)
pangdf <- data.frame(count = as.vector(pangc)/sum(pangc), lineage = names(pangc), df = "pang")
pangdf <- pangdf[pangdf$count > 0.005, ]

allc <- table(alldata$pangolin_lineage)
alldf <- data.frame(count = as.vector(allc)/sum(allc), lineage = names(allc), df = "all")
alldf <- alldf[alldf$lineage %in% pangdf$lineage, ]

both <- rbind(pangdf, alldf)

ggplot(both, aes(x = lineage, y = count, colour = df)) + 
    geom_point() +
    labs(title = "Samples versus all lineages")
