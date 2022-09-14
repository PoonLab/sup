# Download files based on files marked in the NCBI SRA export run info dialog
library(here)

# Load current info (checks to see which files are already done)
acc_file <- here("data", "SraRunInfo_updated.csv")
if(!file.exists(acc_file)) {
    newacc <- read.csv(here("data", "SraRunInfo.csv"))
} else {
    newacc  <- read.csv(acc_file)
}

# If Status column doesn't exist, then assume none have been downloaded yet
if(!"Status" %in% names(newacc)) newacc$Status <- "Not Run"

if(!dir.exists(here("data", "sam"))) dir.create(here("data", "sam"))
if(!dir.exists(here("data", "badsam"))) dir.create(here("data", "badsam"))

# Because export PATH doesn't seem to work well from system call
samdump <- "sam-dump"

# Loop through rows in the csv, do a bunch of error checking
for (i in (1:nrow(newacc))){
    cat("\n")
    t0 <- Sys.time()
    print(t0)
    
    if (newacc$Status[i] != "Not Run") {
        print(paste0("Row ", i, " was already marked ", newacc$Status[i], 
            ", moving on to next row."))
        next
    }
    
    thisnom <- newacc$Run[i]
    thisfil <- here("data", "sam", paste0(thisnom, ".sam"))
    thisbad <- here("data", "badsam", paste0(thisnom, ".sam"))
    
    # Avoid downloading again
    if(!file.exists(thisfil)){
        com <- paste0(samdump, " ", thisnom, " > ", thisfil)
        print(com)
        system(com)
        print(paste0('Download took ', 
            round(difftime(Sys.time(), t0, units = "mins"), 2),
            " minutes."))
    } 
    
    # Check again - if still no file, download must have failed.
    if (!file.exists(thisfil)) {
        print("Download Failed")
        newacc$Status[i] <- "Download Failed"
        write.csv(newacc, acc_file)
        next
    }
    
    # Test the first 1000 lines, rather than reading in the whole file
    firstfew <- readLines(thisfil, 1000)
    if (!(length(firstfew) == 1000)) {
        newacc$Status[i] <- "Incomplete File"
        print("Incomplete File")
        file.rename(thisfil, thisbad)
        write.csv(newacc, acc_file)
        next
    }
    
    # Check if the thousandth line is NA
    testline <- firstfew[1000]
    if(is.na(testline)){
        newacc$Status[i] <- "NA Line"
        print("NA Line")
        file.rename(thisfil, thisbad)
        write.csv(newacc, file = acc_file)
        next
    }

    # If it's a valid file, check the next lines
    if(substr(testline, 1, 1) == "@"){
        nextfew <- readLines(thisfil, 2000)
        testline <- nextfew[2000]
    }
    
    # Check whether cigar string exists
    testcigar <- strsplit(testline, split = "\t")[[1]][6]
    if(length(testcigar) < 1 | is.na(testcigar)) {
        newacc$Status[i] <- "No Cigar"
        print("No Cigar")
        file.rename(thisfil, thisbad)
        write.csv(newacc, file = acc_file)
        next
    }

    # Check whether the cigar string is more than just "*"
    if (nchar(testcigar) < 3) {
        newacc$Status[i] <- "No Cigar"
        print(testcigar)
        print("No Cigar")
        file.rename(thisfil, thisbad)
        write.csv(newacc, file = acc_file)
        next
    }
    
    # If all the error checks have passed, update csv
    newacc$Status[i] <- "Complete"
    write.csv(newacc, file = acc_file)
}




