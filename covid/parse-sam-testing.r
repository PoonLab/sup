source("covid/parse-sam.r")


t0 <- Sys.time()
test1 <- parse.sam("covid/SRR13020990_small.sam")
parsetime <- difftime(Sys.time(), t0, units = "mins")
parsetime

timeperline <- parsetime/500 # this file has 500 lines

# Total time for the file (150,000 lines), in hours
paste(round(as.numeric(timeperline)*147654/60, 3), "hours")
