source("covid/parse-sam.r")


t0 <- Sys.time()
test1 <- parse.sam("covid/SRR13020990_small.sam", time.step = 5)
parsetime <- difftime(Sys.time(), t0, units = "mins")

timeperline <- parsetime/500

# Total time for the file, in hours
paste(round(as.numeric(timeperline)*147654/60, 3), " hours")
