# Download files based on files marked in the NCBI SRA export run info dialog

# Load current info (checks to see which files are already done)
newacc <- read.csv("SraRunInfo_updated.csv")
newsam <- paste0(newacc$Run, ".sam")

# If Status column doesn't exist, then assume none have been downloaded yet
if(!"Status" %in% names(newacc)) newacc$Status <- "Not Run"

# Find all files that have already been downloaded (just to be sure)
alldl <- c(list.files(pattern = "*\\.sam"), list.files("paired", pattern = "*\\.sam"), 
    list.files("unpaired", pattern = "*\\.sam"), list.files("badsam", pattern = "*\\.sam"))


# Because export PATH doesn't seem to work well from system call
samdump <- "~/Downloads/sratoolkit.2.10.8-ubuntu64/bin/sam-dump"

# Loop through rows in the csv, do a bunch of error checking
for (i in (1:length(newsam))){
    cat("\n")
    t0 <- Sys.time()
    print(t0)
    
    if (newacc$Status[i] != "Not Run") {
        print(paste0("Row ", i, " was already marked ", newacc$Status[i], 
            ", moving on to next row."))
        next
    }
    
    thisnom <- newacc$Run[i]
    # Folders are named paired and unpaired
    thispir <- ifelse(newacc$LibraryLayout[i] == "PAIRED", "paired", "unpaired")
    thisfil <- paste0(thispir, "/", thisnom, ".sam")
    thisfas <- paste0(thispir, "/", thisnom, ".fasta")
    
    # Avoid downloading again
    if(!file.exists(thisfil)){
        com <- paste0(samdump, " ", thisnom, " > ", thisfil)
        print(com)
        system(com)
        print(paste0('Download took ', 
            round(difftime(Sys.time(), t0, units = "mins"), 2),
            " minutes."))
    } 
    
    if (!file.exists(thisfil)) {
        print("Download Failed")
        newacc$Status[i] <- "Download Failed"
        write.csv(newacc, file = "SraRunInfo_updated.csv")
        next
    }
    
    # Test the first 1000 lines, rather than reading in the whole file
    firstfew <- readLines(thisfil, 1000)
    if (!(length(firstfew) == 1000)) {
        newacc$Status[i] <- "Incomplete File"
        print("Incomplete File")
        write.csv(newacc, file = "SraRunInfo_updated.csv")
        next
    }
    
    # Check if the thousandth line is NA
    testline <- firstfew[1000]
    if(is.na(testline)){
        newacc$Status[i] <- "NA Line"
        print("NA Line")
        write.csv(newacc, file = "SraRunInfo_updated.csv")
        next
    }
    # If it's a valid file, check the next lines
    if(substr(testline, 1, 1) == "@"){
        nextfew <- readLines(thisfil, 5000)
        testline <- nextfew[5000]
    }
    
    testcigar <- strsplit(testline, split = "\t")[[1]][6]
    if(length(testcigar) < 1 | is.na(testcigar)) {
        newacc$Status[i] <- "No Cigar"
        print("No Cigar")
        write.csv(newacc, file = "SraRunInfo_updated.csv")
        next
    }
    if (nchar(testcigar) < 3) {
        newacc$Status[i] <- "No Cigar"
        print(testcigar)
        print("No Cigar")
        write.csv(newacc, file = "SraRunInfo_updated.csv")
        next
    }
    
    # If all the error checks have passed, download and update csv
    system(paste0("samtools fasta ", thisfil, " > ", thisfas))
    newacc$Status[i] <- "Complete"
    write.csv(newacc, file = "SraRunInfo_updated.csv")
}




