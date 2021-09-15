# Download Files

dls <- read.csv("SRA_downloads.csv", 
    stringsAsFactors = FALSE)

for (i in seq_len(nrow(dls))) { 
    cat("\n")
    t0 <- Sys.time()
    print(t0)

    thisnom <- dls$sra[i]
    thisfil <- paste0(thisnom, ".sam")
    thisfas <- paste0(thisnom, ".fasta")

    if (!dls$status[i] %in% c("Not Run")) {
        print(paste0("Row ", i, " was already marked ", dls$status[i], 
            ", moving on to next row."))

        if (file.exists(thisfil)) {
            file.remove(thisfil)
        }
        next
    }


    
    # Avoid downloading again
    if(!file.exists(thisfil)){
        com <- paste0("sam-dump ", thisnom, " > ", thisfil)
        print(com)
        system(com)
        print(paste0('Download took ', 
            round(difftime(Sys.time(), t0, units = "mins"), 2),
            " minutes."))
    } 


    if (!file.exists(thisfil)) {
        print("Download Failed")

        dls <- read.csv("SRA_downloads.csv", 
            stringsAsFactors = FALSE)
        dls$status[i] <- "Download Failed"
        write.csv(dls, file = "SRA_downloads.csv", row.names = FALSE)
        next
    }
    
    # Test the first 1000 lines, rather than reading in the whole file
    firstfew <- readLines(thisfil, 1000)
    if (!(length(firstfew) == 1000)) {
        print("Incomplete File")

        dls <- read.csv("SRA_downloads.csv", 
            stringsAsFactors = FALSE)
        dls$status[i] <- "Incomplete File"
        write.csv(dls, file = "SRA_downloads.csv", row.names = FALSE)
        if (file.exists(thisfil)) {
            file.remove(thisfil)
        }
        next
    }
    
    # Check if the thousandth line is NA
    testline <- firstfew[1000]
    if(is.na(testline)){
        print("NA Line")

        dls <- read.csv("SRA_downloads.csv", 
            stringsAsFactors = FALSE)
        dls$status[i] <- "NA Line"
        write.csv(dls, file = "SRA_downloads.csv", row.names = FALSE)
        if (file.exists(thisfil)) {
            file.remove(thisfil)
        }
        next
    }
    # If it's a valid file, check the next lines
    if(substr(testline, 1, 1) == "@"){
        nextfew <- readLines(thisfil, 5000)
        testline <- nextfew[5000]
    }
    
    testcigar <- strsplit(testline, split = "\t")[[1]][6]
    if(length(testcigar) < 1 | is.na(testcigar)) {
        print("No Cigar")

        dls <- read.csv("SRA_downloads.csv", 
            stringsAsFactors = FALSE)
        dls$status[i] <- "No Cigar"
        write.csv(dls, file = "SRA_downloads.csv", row.names = FALSE)
        if (file.exists(thisfil)) {
            file.remove(thisfil)
        }
        next
    }
    if (nchar(testcigar) < 3) {
        print(testcigar)
        print("No Cigar")

        dls <- read.csv("SRA_downloads.csv", 
            stringsAsFactors = FALSE)
        dls$status[i] <- "No Cigar"
        write.csv(dls, file = "SRA_downloads.csv", row.names = FALSE)
        if (file.exists(thisfil)) {
            file.remove(thisfil)
        }
        next
    }
    
    # If all the error checks have passed, update csv
    dls <- read.csv("SRA_downloads.csv", 
        stringsAsFactors = FALSE)
    dls$status[i] <- "Complete"
    write.csv(dls, file = "SRA_downloads.csv", row.names = FALSE)

    command <- paste0("../src/parse-sam-c/sam2aln ", thisfil, " 3")
    print(command)
    system(command)
    #file.remove(thisfil)
    print("File was complete; uncertainty matrix generated; file deleted.")
}





