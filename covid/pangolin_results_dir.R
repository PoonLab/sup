# Packages that Art hates
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

# Read in CSV files
csvs <- list.files("data/pangolineages/", pattern = "*.csv")
# Remove any copies
csvs <- csvs[!grepl("-1", csvs)]
# Proper names
csvs <- paste0("data/pangolineages/", csvs)

# Bring them into one data frame
lins <- bind_rows(lapply(csvs, read.csv))

# Taxon is encoded as _ACCSESSIONNUMBER.ID, split into ACCESSIONNUMBER and ID
lins <- lins %>% 
    separate(col = "taxon", sep = "\\.", into = c("taxon", "sample")) %>% 
    mutate(taxon = str_replace(taxon, "\\_", ""))

#### Visualize the uncertainty in the base calls ----
taxons <- unique(lins$taxon)
length(taxons)

# Bar plot (Pareto) for which lineages were called
par(mfrow = c(3, 4))
for(i in 1:length(unique(lins$taxon))){
    thistab <- table(lins$lineage[lins$taxon == taxons[i]])
    barplot(sort(thistab, decreasing = TRUE), las = 2)
}



#### Evaluate pangolin's bootstrap confidence ----
# Find out the modal lineage, look at all of the bootstrap probability,
# then find out how many actually were that lineage

boots <- vector(mode = "list", length = length(taxons))
par(mfrow = c(3,4))
for(i in 1:length(taxons)){
    # Find modal lineage
    thistab <- table(lins$lineage[lins$taxon == taxons[i]])
    maxtab <- names(thistab)[which.max(thistab)]
    # Record "probability" for all sequences labelled with this lineage
    boots[[i]] <- lins$probability[lins$taxon == taxons[i] & 
            lins$lineage == maxtab]
    
    # histogram of probabilities of that lineage
    hist(boots[[i]], breaks = seq(0, 1, 0.05), xlim = c(0, 1))
    # Add a line for the number of genomes that actually had that lineage
    abline(v = mean(lins$lineage[lins$taxon == taxons[i]] == maxtab),
        col = 2, lwd = 3)
}
# Conclusions:
#   The "probability" does NOT account for sequence uncertainty
#   Many sequences had 100% probability for a given lineage,
#       but this is entirely unconnected to the actual probability

# The plots above would be better if I had the consensus sequence, rather
# than modal lineage. 

