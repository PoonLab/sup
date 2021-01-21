library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

csvs <- list.files("data/pangolineages/")
csvs <- csvs[!grepl("-1", csvs)]
csvs <- paste0("data/pangolineages/", csvs)

lins <- bind_rows(lapply(csvs, read.csv))

lins <- lins %>% 
    separate(col = "taxon", sep = "\\.", into = c("taxon", "sample")) %>% 
    mutate(taxon = str_replace(taxon, "\\_", ""))
sort(table(lins$lineage))

taxons <- unique(lins$taxon)
length(taxons)

par(mfrow = c(3, 4))
for(i in 1:length(unique(lins$taxon))){
    thistab <- table(lins$lineage[lins$taxon == taxons[i]])
    #thistab <- thistab[thistab > 1]
    barplot(sort(thistab, decreasing = TRUE), las = 2)
}



#### Evaluate pangolins bootstrap confidence ----
boots <- vector(mode = "list", length = length(taxons))
par(mfrow = c(3,4))
for(i in 1:length(taxons)){
    thistab <- table(lins$lineage[lins$taxon == taxons[i]])
    maxtab <- names(thistab)[which.max(thistab)]
    boots[[i]] <- lins$probability[lins$taxon == taxons[i] & 
            lins$lineage == maxtab]
    
    hist(boots[[i]], breaks = seq(0, 1, 0.05), xlim = c(0, 1))
    abline(v = mean(lins$lineage[lins$taxon == taxons[i]] == maxtab),
        col = 2, lwd = 3)
}

