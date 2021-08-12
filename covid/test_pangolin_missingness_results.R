# Absolute path because the files are too large for github
path <- "/mnt/BCC20BCCC20B8A3A/Fasta Files/sra-downloads/files/pangolin_results/"
finished <- list.files(path)
#length(finished)

par(mfrow = c(4,4))
for(i in seq_along(finished)){
    pang <- read.csv(paste0(path, finished[i]))
    pang$taxon <- gsub("\\_", "", pang$taxon)
    pang$missing <- as.numeric(sapply(strsplit(pang$taxon, "\\."), 
        function(x) x[1]))

    boxplot(pang$probability ~ pang$missing, main = finished[i],
        las = 2, xlab = NULL, ylim = c(0, 1))
}

# barplot(table(pang$probability, pang$missing))


# TODO
    # For ERR4693034 and ERR4692364, investigate more values
        # Fewer replicates, more proportions
    # Try adding a logistic regression curve
    # Investigate: Maximum for SRR* files is less than 1? Why?


# The new ones:
path <- "/mnt/BCC20BCCC20B8A3A/Fasta Files/sra-downloads/files/pangolin_results/"
finished <- list.files(path, pattern = "*-2.csv")
#length(finished)

par(mfrow = c(1,2))
for(i in seq_along(finished)){
    pang <- read.csv(paste0(path, finished[i]))
    pang$taxon <- gsub("\\_", "", pang$taxon)
    pang$missing <- as.numeric(sapply(strsplit(pang$taxon, "\\."), 
        function(x) x[1]))
    
    #boxplot(pang$probability ~ pang$missing, main = finished[i],
    #    las = 2, xlab = NULL, ylim = c(0, 1))
    plot(pang$probability ~ pang$missing, main = finished[i])
    
    plm <- glm(I(probability > 0.5) ~ missing, family = binomial, data = pang)
    xseq <- seq(min(pang$missing), max(pang$missing), len = 100)
    lines(x=xseq, y = predict(plm, newdata = list(missing = xseq), type = "response"))
}

