source("covid/parse-sam.r")

t0 <- Sys.time()
test1 <- parse.sam("data/SRR13020990_small.sam")
parsetime <- difftime(Sys.time(), t0, units = "mins")
parsetime

timeperline <- parsetime/500 # this file has 500 lines

# Total time for the file (150,000 lines), in hours
paste(round(as.numeric(timeperline)*147654/60, 3), "hours")

source("covid/parse-sam_dev.r")
tn <- Sys.time()
test2 <- parse.sam_deprecated("data/SRR13020990_small.sam")
parsetime <- difftime(Sys.time(), tn, units = "mins")
parsetime

timeperline <- parsetime/500 # this file has 500 lines

# Total time for the file (150,000 lines), in hours
paste(round(as.numeric(timeperline)*147654/60, 3), "hours")

all.equal(test1, test2)

tn <- Sys.time()
test2 <- parse.sam.dt("data/SRR13020990_small.sam")
parsetime <- difftime(Sys.time(), tn, units = "mins")
parsetime

timeperline <- parsetime/500 # this file has 500 lines





# Paired read testing
test3 <- parse.sam("data/SRR13592146_small.sam")

