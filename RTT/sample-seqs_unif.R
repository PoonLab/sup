# Sample uniformly from weeks
# Libraries and Constants -----------------------
start_time <- Sys.time()
library(dplyr)

samples_per_week <- 1
args <- commandArgs(trailing = TRUE)
if("-N" %in% args){
	samples_per_week <- as.numeric(args[which(args == "-N") + 1])
}

# Load and Filter Data --------------------------
desc <- read.csv("sequences_descr_wk.csv")
# Calculate number of lines per genome (each row is 60 nucleotides)
desc <- desc %>%
	arrange(RowNum) %>%
	mutate(Lines = (Length %/% 60) + 1)

week_counts <- desc %>% 
	group_by(yr, wk) %>%
	count()

desc <- left_join(desc, week_counts, by = c("yr", "wk")) %>% 
	filter(n >= samples_per_week)


# Gather Information ----------------------------
set.seed(2112); seq_per_week <- desc %>% 
	filter(paste(yr, wk, sep = "-") != "2019-18") %>%
	group_by(yr, wk) %>%
	mutate(
		sampled = RowNum %in% sample(RowNum, 
			size = min(n(), samples_per_week),
			replace = FALSE)
	) %>% 
	filter(sampled) %>%
	ungroup() %>% 
	arrange(RowNum)
seq_per_week

if (FALSE) { # check that it's uniform
	library(ggplot2)
	ggplot(seq_per_week) + 
		aes(x = yr + wk/53) + 
		geom_bar()

	hist(seq_per_week$Length)
}




t0 <- Sys.time()
if ("sampled_sequences3.fasta" %in% list.files()) {
	file.remove("sampled_sequences.fasta")
}
file.create("sampled_sequences.fasta")

biggun <- file("sequences.fasta", "r")
for (i in seq_len(nrow(seq_per_week))) {
	if (!i %% 25) {
		print(paste0("i=", i, "/", nrow(seq_per_week), ", ", 
			seq_per_week$Accession[i]))
	}
	start <- seq_per_week$RowNum[i]

	if (i == 1) {
		throw_out <- readLines(biggun, n = start - 1)
	} else {
		so_far <- seq_per_week$RowNum[i - 1] + seq_per_week$Lines[i-1]
		n_throw <- seq_per_week$RowNum[i] - so_far
		# Fudging n_throw so we undershoot the line number
		# Several of the seq_per_week$Lines are off by 3 or 4
		n_throw <- n_throw - 10
		if(n_throw > 0) {
			throw_out <- readLines(biggun, n = n_throw)
		}
	}

	found_it <- FALSE
	while (!found_it) {
		lines1 <- readLines(biggun, n = 1)
		counter <- 0
		if (grepl(pattern = ">", x = lines1)) {
			counter <- counter + 1
			if(counter >= 50) {
				print("Too many reads")
				break
			}

			acc <- strsplit(lines1, "\\|")[[1]][1]
			acc <- gsub(pattern = ">", replacement = "", x = acc)
			#acc <- trimws(acc)
			#print(acc)

			if (acc != seq_per_week$Accession[i]) {
				print("Uh oh")
				print("Accession number mismatch.")
				print(paste0("acc:", acc))
				print(paste0("Accession[i - 1]:", 
					seq_per_week$Accession[i - 1]))
				print(paste0("Accession[i]:", 
					seq_per_week$Accession[i]))
				print(paste0("Accession[i + 1]:", 
					seq_per_week$Accession[i + 1]))
				break
			} else {
				found_it <- TRUE
			}

			write(lines1, 
				file = "sampled_sequences.fasta", 
				append = TRUE)
			lines <- readLines(biggun, n = seq_per_week$Lines[i])
			write(lines, 
				file = "sampled_sequences.fasta", 
				append = TRUE)
		}	
	}
}
close(biggun)
print(Sys.time() - t0)

cat("\nTotal time:", 
	difftime(Sys.time(), start_time, units = "mins"), 
	"mins \n")



if (FALSE) { # Deprecated (slower way)
	# Create sampled fasta file ---------------------
	# Loop through sampled values, append to new file
	# sed -n '1001,1500p;1500q' sequences.fasta >> sampled_sequences.fasta

	if ("sampled_sequences.fasta" %in% list.files()) {
		file.remove("sampled_sequences.fasta")
	}
	file.create("sampled_sequences.fasta")

	for (i in seq_len(nrow(seq_per_week))) {
		start <- seq_per_week$RowNum[i]
		end <- start + seq_per_week$Lines[i]
		quitline <- end + 1

		sed <- paste0("sed -n '",
			start, ",",
			end, "p;",
			quitline, "q'",
			" sequences.fasta >> sampled_sequences.fasta")
		cat("\n", i, " | ", sed)
		# Run it!
		system(sed)
	}
}
