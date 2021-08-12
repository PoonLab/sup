# Quality score by sample date
library(lubridate)

sams <- list.files("samDONE")
sams <- sams[!grepl(x = sams, pattern = "insertions")]
sams <- sapply(sams, function(x) strsplit(x, split = "\\.")[[1]][1])
sams <- sapply(sams, function(x) strsplit(x, split = "-")[[1]][2])
sams


meta <- read.csv("sequences_descr_mt.csv")
meta <- meta[!is.na(meta$SRAAccession), ]

any(!sams %in% meta$SRAAccession)

s1 <- sams[1]

date1 <- meta[meta$SRAAccession == s1, "date"]

sam1 <- readRDS(paste0("samDONE/S-", s1, ".RDS"))
sam1 <- sam1[2:5]

datequal <- data.frame(acc = sams,
    coverage = NA, date = NA, len = NA, summax = NA)

for (i in seq_len(nrow(datequal))) {
    sami <- readRDS(paste0("samDONE/S-", sams[i], ".RDS"))
    datequal$len[i] = nrow(sami)
    datequal$summax[i] <- sum(sapply(seq_along(nrow(sami)),
        function(x) max(sami[x, ])))
    datequal$coverage[i] <- sum(sami)
    datequal$date[i] <- meta[meta$SRAAccession == sams[i], "date"]
}
datequal$date <- ymd(datequal$date)

plot(coverage ~ date, data = datequal); abline(lm(coverage ~ date, data = datequal))

plot(summax ~ date, data = datequal); abline(lm(summax ~ date, data = datequal))


