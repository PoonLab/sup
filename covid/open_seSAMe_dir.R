# To run in background (esp. on rei): 
# nohup Rscript open_seSAMe_dir [dirname] &

source("parse-sam.r")

args <- commandArgs(TRUE)

# This file assumes a particular directory structure.
# SAM files must either be in current_dir or current_dir/samTODO
samfile <- list.files(path = args[1], pattern = "*\\.sam")
if(length(samfile) == 0){
    stop("no files in specified directoy")
}

if(!dir.exists("samDONE")) dir.create("samDONE")

print(length(samfile))

for(f in seq_along(samfile)){
    cat("\n")
    t0 <- Sys.time()
    thisfilename <- paste0(args[1], "/", samfile[f])
    print(thisfilename)
    
    thisacs <- substr(samfile[f], 1, nchar(samfile[f]) - 4)
    thisRDS <- paste0("samDONE/S-", thisacs, ".RDS")
    
    if(file.exists(thisRDS)) next
    
    try({ # don't break loop if there's an error
        test1 <- parse.sam.dt(thisfilename, nc = 1)
        
        saveRDS(test1, file = thisRDS)
        
        tt <- round(difftime(Sys.time(), t0, units = "hours"), 3)
        msg1 <- paste0(thisacs, " completed in ", tt, " hours.")
        msg2 <- paste0("Result has ", ncol(test1), " cols and ", nrow(test1), " rows.")
        msg <- paste(msg1, msg2, collapse = "/n")
        print(msg)
        writeLines(msg, paste0("open_seSAMe_info-", thisacs, ".txt"))
    })
}

print("Done!")

