# Data Directory



`downloader.R` downloads and validates all sequences used in this study, as described by the `SraRunInfo.csv` file. Depends on the `sam-dump` binary from `sratoolkit` and uses the the R package called `here` so that it can be run from any directory in the repo.

As a by-product, it creates the file `SraRunInfo_updated.csv`, which has a column that reports the result of the 

Any file with errors - incomplete files, missing CIGAR strings, etc. - are sorted into the `badsam` directory. All successful downloads are dumped into the `sam` directory. The script makes these directories if they do not exist.

As a by-product, it creates the file `SraRunInfo_updated.csv`, which has a column that reports any failures during the downloading process.

To run this on other SRA runs, simply changing the SraRunInfo file is sufficient. If you already have SAM files, then you should skip this step.

---

`run_all.sh` runs `parse-sam-c` on all sam files in the `sam` directory. This script must be run from the root directory.

This results in a folder `output` which contains the output of this function.

I started it at approximately 2:30pm on September 15th, we'll see how long it takes before the first file is completed. 






