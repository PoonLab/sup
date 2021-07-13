args <- commandArgs(TRUE)
infile <- 
seq_names <- readLines(args[1])
# Remove the ">" at the start of every line
seq_names <- substr(seq_names, 2, 1000)
seq_dates <- sapply(strsplit(seq_names, split = "\\|"), function(x) x[3])
write.csv(data.frame(name = seq_names, date = seq_dates), 
	file = args[2], row.names = FALSE, quote = FALSE)
