files <- list.files(pattern = '(ERR|SRR)')
files <- files[grepl(files, pattern = "csv")]

for(f in files) {
	print(f)
	saveRDS(read.csv(f), file = paste0('samDONE/S-', sub('csv', 'RDS', f)))
	if ("-d" %in% commandArgs()) file.remove(f)
}
