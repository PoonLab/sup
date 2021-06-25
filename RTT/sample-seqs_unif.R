# Sample uniformly from weeks
# Libraries and Constants -----------------------
suppressPackageStartupMessages({
	library(dplyr)
	library(lubridate)
})

samples_per_month <- 3
args <- commandArgs(trailing = TRUE)
if("-N" %in% args){
	samples_per_month <- as.numeric(args[which(args == "-N") + 1])
}

# Load and Filter Data --------------------------
desc <- read.csv("sequences_descr_mt.csv")
# Calculate number of lines per genome (each row is 60 nucleotides)
desc <- desc %>%
	arrange(RowNum) %>%
	mutate(Lines = (Length %/% 60) + 1,
		mt = month(ymd(date)))

week_counts <- desc %>% 
	group_by(yr, mt) %>%
	count()

desc <- left_join(desc, week_counts, by = c("yr", "mt")) %>% 
	filter(n >= samples_per_month) %>%
	filter(nchar(desc$SRAAccession) > 0)


# Gather Information ----------------------------
set.seed(2); seq_per_week <- desc %>% 
	filter(paste(yr, wk, sep = "-") != "2019-18") %>%
	group_by(yr, mt) %>%
	mutate(
		sampled = RowNum %in% sample(RowNum, 
			size = min(n(), samples_per_month),
			replace = FALSE)
	) %>% 
	filter(sampled) %>%
	ungroup() %>% 
	arrange(RowNum)
#seq_per_week

if (FALSE) { # check that it's uniform
	library(ggplot2)
	ggplot(seq_per_week) + 
		aes(x = yr + mt/13) + 
		geom_bar()

	# Sanity check
	hist(seq_per_week$Length)
}

defs <- seq_per_week$def
defs <- sapply(defs, function(x) strsplit(x, split = ":>")[[1]][2])

writeLines(defs, con = "sampled_seqs.txt")


# SRA Accession Numbers

sra <- seq_per_week$SRAAccession

writeLines(sra, con = "sampled_SRA.txt")

