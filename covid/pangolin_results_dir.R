# Packages that Art hates
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

args <- commandArgs(TRUE)
dirich <- "-d" %in% args

# Read in CSV files
csvs <- list.files("data/pangolineages/", 
    pattern = ifelse(dirich, "*_d.csv", "*.csv"))
# Remove any copies
csvs <- csvs[!grepl("-1", csvs)]
# Proper names
csvs <- paste0("data/pangolineages/", csvs)

# Bring them into one data frame
lins <- bind_rows(lapply(csvs, read.csv))

# Taxon is encoded as _ACCSESSIONNUMBER.ID, split into ACCESSIONNUMBER and ID
lins <- lins %>% 
    separate(col = "taxon", sep = "\\.", 
        into = c("taxon", "sample")) %>% 
    mutate(taxon = str_replace(taxon, "\\_", ""))

#### Visualize the uncertainty in the base calls ----
taxons <- unique(lins$taxon)
length(taxons)

# Bar plot (Pareto) for which lineages were called
png(ifelse(dirich, "figures/pareto_d.png", "figures/pareto.png"), 
    height = 1000, width = 1000, units = "px")
par(mfrow = c(6, 6))
for(i in 1:length(unique(lins$taxon))){
    thistab <- table(lins$lineage[lins$taxon == taxons[i]])
    called <- lins$lineage[lins$taxon == taxons[i] & 
        lins$sample == 0]
    print(called)
    
    if(length(called) == 1) {
        colours <- rep(1, length(thistab))
        colours[names(thistab) == called] <- 2
    } else {
        colours <- "lightgrey"
    }
    
    barplot(sort(thistab, decreasing = TRUE), 
        las = 2, col = colours, border = NA)
}
dev.off()

# Info for the abstract
summs <- lins %>% 
    group_by(taxon) %>%
    summarise(maxperc = mean(lineage == names(sort(table(lineage), 
        decreasing = TRUE))[1]),
        uniques = length(unique(lineage)),
        minpango = min(probability),
        maxpango = max(probability),
        menpango = mean(probability),
        max = names(sort(table(lineage), decreasing = TRUE))[1])

print("summary info")
print(summs)
print(1 - mean(summs$maxperc); 1 - mean(summs$menpango))

#### Evaluate pangolin's bootstrap confidence ----
# Find out the modal lineage, look at all of the bootstrap probability,
# then find out how many actually were that lineage

png(ifelse(dirich, "figures/bar_d.png", "figures/bar.png"), 
    height = 1000, width = 1000, units = "px")
boots <- vector(mode = "list", length = length(taxons))
par(mfrow = c(6,6))
for(i in 1:length(taxons)){
    # Find modal lineage
    thistab <- table(lins$lineage[lins$taxon == taxons[i]])
    maxtab <- names(thistab)[which.max(thistab)]
    # Record "probability" for all sequences labelled with this lineage
    boots[[i]] <- lins$probability[lins$taxon == taxons[i] & 
            lins$lineage == maxtab]
    
    # histogram of probabilities of that lineage
    hist(boots[[i]], breaks = seq(0, 1, 0.05), 
        xlim = c(0, 1))
    # Add a line for the number of genomes that actually had that lineage
    abline(v = mean(lins$lineage[lins$taxon == taxons[i]] == maxtab),
        col = 2, lwd = 3)
}
dev.off()
# Conclusions:
#   The "probability" does NOT account for sequence uncertainty
#   Many sequences had 100% probability for a given lineage,
#       but this is entirely unconnected to the actual probability

# The plots above would be better if I had the consensus sequence, rather
# than modal lineage. 

gglist <- list()
for(i in 1:length(taxons)){
    pang <- lins[lins$taxon == taxons[i], ]
    pang2 <- lapply(unique(pang$lineage), function(x) {
        data.frame(lineage = x, 
            prop = mean(pang$lineage == x))
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
    
    gglist[[i]] <- pang2 %>% 
        # round to nearest 0.5
        #mutate(prop = round(prop*2, 1)/2,
        #    probability = round(probability*2, 1)/2) %>% 
        group_by(prop, probability, lineage) %>% 
        summarise(count = n(), .groups = "drop") %>% 
        ggplot() + theme_bw() + 
        aes(x = prop, y = probability, colour = lineage, 
            label = count) + 
        geom_text() + 
        theme(legend.position = "none") +
        annotate("text", x = pangtab$prop, y = 1, 
            label = pangtab$lineage,
            hjust = 0.5, vjust = -1) +
        labs(x = "Proportion of Lineage",
            y = "Bootstrap Probability",
            title = NULL) +
        scale_x_continuous(breaks = seq(0,1,0.1)) +
        scale_y_continuous(breaks = seq(0,1.1,0.1)) +
        coord_cartesian(ylim = c(0, 1.1)) + 
        geom_abline(slope = 1, intercept = 0)
}
cowplot::plot_grid(plotlist = gglist)
ggsave(ifelse(dirich, "figures/ggscatter_p.png", "figures/ggscatter.png"))












# Temp testing code
pang <- lins[lins$taxon == taxons[1], ]
called <- pang$lineage[pang$sample == 0][1]
pangtab <- sort(table(pang$lineage), decreasing = TRUE)
colvec <- rep("grey", length(pangtab))
colvec[which(names(pangtab) == called)] <- "red"
n <- sum(pangtab > 200)
barlabx <- c(0, cumsum(pangtab[1:(n-1)])) + pangtab[1:n]/2
barplot(as.matrix(pangtab), col = colvec, hori = TRUE)
text(barlabx, 0.7, names(pangtab_small))


#par.old <- par()
png("figures/stacked-lineages.png", 
    width = 4, height = 8, units = "in", res = 300)
par(mfrow = c(15, 1), mar = c(1, 0.5, 1.5, 0.5))
for (i in 1:38) {
    pang <- lins[lins$taxon == taxons[i], ]
    called <- pang$lineage[pang$sample == 0][1]
    pangtab <- sort(table(pang$lineage), decreasing = TRUE)

    colvec <- rep("grey", length(pangtab))
    colvec[which(names(pangtab) == called)] <- "red"

    n <- sum(pangtab > 50)
    if (n > 1) {
        barlabx <- c(0, cumsum(pangtab[1:(n-1)])) + 
            pangtab[1:n]/2

        barplot(as.matrix(pangtab), 
            col = colvec, hori = TRUE, axes = FALSE)
        text(barlabx, 0.7, names(pangtab)[1:n])
        mtext(side = 3, text = taxons[i], cex = 0.75, las = 1)
        pretty_labels = seq(0, sum(pangtab),
                by = ifelse(sum(pangtab) < 2000, 100, 1000))
        mtext(side = 1, 
            at = pretty_labels,
            text = pretty_labels,
            line = 0,
            cex = 0.75
        )
    } else {
        cat("Taxon", taxons[i], "had", length(pangtab),
            "unique calls, with largest accounting for", 
            pangtab[1], "lineage calls, with second place", 
            pangtab[2], ". First place was",
            ifelse(names(pangtab)[1] == called, "", "not"),
            "the conseq call.\n")
    }
}
dev.off()









