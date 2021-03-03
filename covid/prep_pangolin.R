# # Commands for shell (can't be run from RStudio because RStudio uses the wrong path)
# conda activate pangolin
# pangolin data/SRR13020989_sampled.fasta --outfile data/SRR13020989_pangolineages.csv
# # for 1000 sequences, took 30 seconds

done <- list.files("data/pangolineages", pattern = "*\\.csv")
done <- sapply(strsplit(done, split = "_"), function(x) x[1])

ready <- list.files("data/sampled_covid", pattern = "*\\.fasta")
ready <- sapply(strsplit(ready, split = "_"), function(x) x[1])

todo <- ready[!ready %in% done]
todo

for(i in 1:length(todo)){
    cat(paste0(
        "pangolin data/sampled_covid/", todo[i], "_sampled.fasta",
        " --outfile data/pangolineages/", todo[i], "_pangolineages.csv"
    ))
    cat("\n")
}; cat('spd-say "Hey devan your code is done"\n\n')

