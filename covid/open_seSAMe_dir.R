# To run in background (esp. on rei): 
# nohup Rscript open_seSAMe_multi.R ERR4305816.sam ERR4085809.sam ERR4204823.sam ERR4307842.sam

source("parse-sam.r")

# This file assumes a particular directory structure.
# SAM files must either be in current_dir or current_dir/samTODO
samfile <- list.files(path = "samTODO", pattern = "*.sam")
if(length(samfile) == 0){
    stop("no files specifed")
}
if(!all(grepl(".sam", samfile))){
    stop("at least one bad file name")
}
# Checks current_dir AND samTODO
if(!all(samfile %in% c(list.files(), list.files("samTODO")))){
    stop("files must be in directory called samTODO")
}

if(!dir.exists("samDONE")) dir.create("samDONE")

print(length(samfile))

for(f in seq_along(samfile)){
    t0 <- Sys.time()
    thisfilename <- samfile[f]
    print(thisfilename)
    
    thisacs <- substr(thisfilename, 1, nchar(thisfilename) - 4)
    
    # Check whether it's current dir or sam TODO
    if(!thisfilename %in% list.files()){
        thisfilename <- paste0("samTODO/", thisfilename, collapse = "")
    }
    
    # Estimate length for timing purposes - first 4 lines are metadata
    # This is approximate (and wrong, but close enough)
    est.length <- length(count.fields(thisfilename)) - 4
    print(est.length)
    
    try({ # don't break loop if there's an error
        test1 <- parse.sam(thisfilename, est.length = est.length, time.step = 2000, paired=FALSE)
        
        saveRDS(test1, file = paste0("samDONE/open_seSAMe-S-", thisacs, ".RDS"))
        
        tt <- round(difftime(Sys.time(), t0, units = "hours"), 3)
        msg1 <- paste0(thisacs, " completed in ", tt, " hours.")
        msg2 <- paste0("Result has ", ncol(test1), " cols and ", nrow(test1), " rows.")
        msg <- paste(msg1, msg2, collapse = "/n")
        print(msg)
        writeLines(msg, paste0("open_seSAMe_info-", thisacs, ".txt"))
    })
}

print("Done!")

