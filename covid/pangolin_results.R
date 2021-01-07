library(ggplot2)

pang <- read.csv("data/SRR13020989_pangolineages.csv")
alldata <- read.csv("data/metadata.tsv", sep = "\t")

barplot(sort(table(pang$lineage), decreasing = TRUE))

sort(table(pang$lineage), decreasing = TRUE)[1:10]

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
