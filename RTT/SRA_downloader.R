# Download Files

dls <- read.csv("SRA_downloads.csv")

for (i in seq_len(nrow(dls))) { 
    cat("\n")
    t0 <- Sys.time()
    print(t0)

    if (dls$status[i] != "Not Run") {
        print(paste0("Row ", i, " was already marked ", dls$status[i], 
            ", moving on to next row."))
        next
    }


    thisnom <- dls$sra[i]
    # Folders are named paired and unpaired
    thisfil <- paste0(thisnom, ".sam")
    thisfas <- paste0(thisnom, ".fasta")
    
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
        dls$status[i] <- "Download Failed"
        write.csv(dls, file = "SraRunInfo_updated.csv")
        next
    }
    
    # Test the first 1000 lines, rather than reading in the whole file
    firstfew <- readLines(thisfil, 1000)
    if (!(length(firstfew) == 1000)) {
        dls$status[i] <- "Incomplete File"
        print("Incomplete File")
        write.csv(dls, file = "SraRunInfo_updated.csv")
        next
    }
    
    # Check if the thousandth line is NA
    testline <- firstfew[1000]
    if(is.na(testline)){
        dls$status[i] <- "NA Line"
        print("NA Line")
        write.csv(dls, file = "SraRunInfo_updated.csv")
        next
    }
    # If it's a valid file, check the next lines
    if(substr(testline, 1, 1) == "@"){
        nextfew <- readLines(thisfil, 5000)
        testline <- nextfew[5000]
    }
    
    testcigar <- strsplit(testline, split = "\t")[[1]][6]
    if(length(testcigar) < 1 | is.na(testcigar)) {
        dls$status[i] <- "No Cigar"
        print("No Cigar")
        write.csv(dls, file = "SraRunInfo_updated.csv")
        next
    }
    if (nchar(testcigar) < 3) {
        dls$status[i] <- "No Cigar"
        print(testcigar)
        print("No Cigar")
        write.csv(dls, file = "SraRunInfo_updated.csv")
        next
    }
    
    # If all the error checks have passed, update csv
    dls$status[i] <- "Complete"
    write.csv(dls, file = "SRA_downloads.csv")
}





