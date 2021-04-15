# To run in background (esp. on rei): 
# nohup Rscript open_seSAMe2.R 

t0 <- Sys.time()
setwd("/home/dbecker7/run_sam")

# Using a custom version of parse.sam from the sung package
source("parse-sam.r")

test1 <- parse.sam("SRR13020989.sam")
#test2 <- parse.sam.mp("~/PostDoc Fasta Files/sra-downloads/files/SRR13020989_2.sam",
#    n.cores = 6)
# The parallel version (parse.sam.mp) recalculates phred at each position. Yes,
# it's parallel, but with a literal metric ton of calculation overhead.
#str(test1)

saveRDS(test1, file = "open_seSAMe2-S.RDS")

tt <- round(difftime(Sys.time(), t0, units = "hours"), 3)
msg1 <- paste0("Completed in ", tt, " hours.")
msg2 <- paste0("Result has ", ncol(test1), " cols and ", nrow(test1), " rows.")
msg <- paste(msg1, msg2, collapse = "/n")
print(msg)
writeLines(msg, "open_seSAM2_info.txt")




