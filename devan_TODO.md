# Roadmap

-   [ ] Update parse.sam to divide paired reads by 2

-   [ ] Validate code on paired reads

-   [ ] Download many more (paired and unpaired) sam files, maybe 30 or 40?

    -   Note that many will fail/timeout or won't be useable
    -   Find links manually
    -   Create sam-dump commands from Rscript file
    -   Download to airlock drive

-   [ ] Send sam files to Rei

    -   Embarassingly embarassingly parallel: just run the script a few times with different directories.

-   [ ] Use scripts from [sam2conseq](https://github.com/PoonLab/sam2conseq/blob/master/sam2conseq.py) to find consensus sequences for sam files

    -   Add these to the files to be sent into pangolin
    -   Potentially, find the reported sequence (as submitted to NCBI) instead

-   [x] Ensure pangolin is working as expected

    -   Try it out on Rei?

-   [ ] Polish the visualizations of uncertainty

    -   Create new ones?
    -   Some way to measure concordance
    -   Many, many sequences

-   [ ] New folder: covid\_src

    -   Populate with existing source-able scripts, ensure outputs are consistent
    -   Add error-checking steps

-   [ ] Make all remaining scripts source-able, with appropriate outputs

-   [ ] Write pipeline script to source the files

    -   Does not need to include downloading scripts
    -   May be best to separate out parsing scripts as well (even with speedups, 30 sam files will take several days)
