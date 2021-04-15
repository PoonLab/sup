source("covid/parse-sam.r")

inFile <- "data/ERR5082711_small.sam" # unpaired
inFile <- "data/ERR5069871.sam" # Paired
inFile <- "/media/devan/Seagate Basic/PostDoc/files/unpaired/blocku1/ERR4085809.sam"

t0 <- Sys.time()
test1 <- parse.sam(inFile)
parsetime <- difftime(Sys.time(), t0, units = "mins")
parsetime

timeperline <- parsetime/500 # this file has 500 lines

# Total time for the file (150,000 lines), in hours
paste(round(as.numeric(timeperline)*147654/60, 3), "hours")

source("covid/parse-sam_dev.r")
tn <- Sys.time()
test2 <- parse.sam_deprecated(inFile)
parsetime <- difftime(Sys.time(), tn, units = "mins")
parsetime

timeperline <- parsetime/500 # this file has 500 lines

# Total time for the file (150,000 lines), in hours
paste(round(as.numeric(timeperline)*147654/60, 3), "hours")

all.equal(test1, test2)

tn <- Sys.time()
test2 <- parse.sam.dt(inFile, nc = 4)
parsetime <- difftime(Sys.time(), tn, units = "mins")
parsetime

timeperline <- parsetime/500 # this file has 500 lines


# From Gopi's function
test4 <- read.csv("matrix.csv")
test4 <- as.matrix(test4)[,-1]

str(test2)
str(test4)

all.equal(test2, test4)
test4
test4 - test2
