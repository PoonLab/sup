# Roadmap

-   [x] Update parse.sam to divide paired reads by 2

-   [ ] Validate code on paired reads

-   [x] Create txt file with, like, 100 NCBI URL links

-   [x] Parse the html, save relevant info (acc/run number, reads, etc.)

-   [x] Keep track of which ones are done

-   [x] Create/run sam-dump commands using system() and tryCatch()

    -   Alternative: learn bash scripting

-   [x] Put sam files in relevant directories

-   [x] Download many more (paired and unpaired) sam files, maybe 30 or 40?

    -   Note that many will fail/timeout or won't be useable

    -   Find links manually

    -   Add links to a text file, read from that file

    -   Keep track of accession/run numbers

    -   Create sam-dump commands from Rscript file

        -   Possibly even run them from Rscript

    -   Download to airlock drive

-   [x] git clone dsup to Rei

-   [x] Send sam files to Rei

    -   Embarassingly embarassingly parallel: just run the script a few times with different directories.
    -   Re-jig open_seSAMe to take directory as argument
    -   Fix output

-   [x] Use samtools to create fasta from sam

-   [ ] Sample from the uncertainty matrices

-   [ ] Add fasta to the files to be sent into pangolin

-   [ ] Ensure pangolin is working as expected

    -   Ask Art about running it on Rei?

-   [ ] Polish the visualizations of uncertainty

    -   Create new ones?
    -   Some way to measure concordance
    -   Many, many sequences

-   [ ] New folder: covid_src

    -   Populate with existing source-able scripts, ensure outputs are consistent
    -   Add error-checking steps

-   [ ] Make all remaining scripts source-able, with appropriate outputs

-   [ ] Write pipeline script to source the files

    -   Does not need to include downloading scripts
    -   May be best to separate out parsing scripts as well (even with speedups, 30 sam files will take several days)
