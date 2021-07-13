
sams <- list.files("samDONE")
sams <- unique(sapply(strsplit(sams, split = "-|_|\\."), function(x) x[2]))

metadata <- read.csv("sequences_descr_mt.csv")

metadata$dl <- FALSE 
for (i in seq_len(nrow(metadata))) {
	testcase <- metadata$SRAAccession[i]
	testcase <- trimws(strsplit(testcase, split = ",")[[1]])
	if(any(testcase %in% sams)) {
		metadata$dl[i] <- TRUE
	}
}

metadl <- metadata[which(metadata$dl), ]

write.csv(metadl, "sequences_descr_mt_downloaded.csv", row.names = FALSE)

sampled_seq_dl <- metadl$def
sampled_seq_dl <- sapply(strsplit(sampled_seq_dl, split = ":>"), function(x) x[2])

writeLines(sampled_seq_dl, "sampled_seq_dl.txt")
