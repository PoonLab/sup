# Absolute path because the files are too large for github
path <- "/mnt/BCC20BCCC20B8A3A/Fasta Files/sra-downloads/files/pangolin_results/"
finished <- list.files(path)
length(finished)

par(mfrow = c(3,3))
for(i in seq_along(finished)){
    pang <- read.csv(paste0(path, finished[i]))
    pang$taxon <- gsub("\\_", "", pang$taxon)
    pang$missing <- sapply(strsplit(pang$taxon, "\\."), function(x) x[1])
    
    boxplot(pang$probability ~ pang$missing, main = finished[i])
}

barplot(table(pang$probability, pang$missing))
