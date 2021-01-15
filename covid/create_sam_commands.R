# This is just a convenience file to save me from
# manually typing in accession numbers as I find them.

# Copy and paste (or drag and drop) links here.

urls1 <- "https://www.ncbi.nlm.nih.gov/sra/SRX9537817[accn]
https://www.ncbi.nlm.nih.gov/sra/SRX9537816[accn]
https://www.ncbi.nlm.nih.gov/sra/ERX4615584[accn]
https://www.ncbi.nlm.nih.gov/sra/ERX4615465[accn]
https://www.ncbi.nlm.nih.gov/sra/ERX4614837[accn]
https://www.ncbi.nlm.nih.gov/sra/ERX4614688[accn]
https://www.ncbi.nlm.nih.gov/sra/ERX4614620[accn]
https://www.ncbi.nlm.nih.gov/sra/ERX4614165[accn]
https://www.ncbi.nlm.nih.gov/sra/ERX4614147[accn]
https://www.ncbi.nlm.nih.gov/sra/ERX4614119[accn]
https://www.ncbi.nlm.nih.gov/sra/SRX9233401[accn]
https://www.ncbi.nlm.nih.gov/sra/ERX4613447[accn]
https://www.ncbi.nlm.nih.gov/sra/ERX4590348[accn]
https://www.ncbi.nlm.nih.gov/sra/ERX4582858[accn]
https://www.ncbi.nlm.nih.gov/sra/ERX4307528[accn]
https://www.ncbi.nlm.nih.gov/sra/ERX4307447[accn]
"


urls <- strsplit(urls1, split = "\n")[[1]]

# Process the URLs
commands <- c()
for(i in seq_along(urls)[-(1:7)]){
    test1 <- xml2::read_html(urls[i])
    test2 <- as.character(test1)
    # The string "run=" only appears at the accession number
    pos <- gregexpr("run=", test2)[[1]]
    # Sanity check
    print(substr(test2, pos + 3, pos + 15))
    # accession numbers are 9-11 characters long. Requires post-processing
    asc <- substr(test2, pos + 4, pos + 14)
    com <- paste0("sam-dump ", asc, " > ", asc, ".sam")
    commands <- c(commands, com)
}
cat(commands, sep = "\n")

# Copying and pasting the output here to modify codes if necessary
" 
sam-dump SRR13092001 > SRR13092001.sam # timed out
sam-dump SRR13092002 > SRR13092002.sam
sam-dump ERR4694498 > ERR4694498.sam
sam-dump ERR4694380 > ERR4694380.sam
sam-dump ERR4758772 > ERR4758772.sam
sam-dump ERR4693605 > ERR4693605.sam
sam-dump ERR4693537 > ERR4693537.sam
sam-dump ERR4693079 > ERR4693079.sam
sam-dump ERR4693061 > ERR4693061.sam
sam-dump ERR4693034 > ERR4693034.sam
sam-dump SRR12762573 > SRR12762573.sam
sam-dump ERR4692364 > ERR4692364.sam
sam-dump ERR4667618 > ERR4667618.sam
sam-dump ERR4664555 > ERR4664555.sam
sam-dump ERR4363387 > ERR4363387.sam
sam-dump ERR4364007 > ERR4364007.sam
"