# Roadmap

-   [ ] Download many more (UNPAIRED) sam files, maybe 30 or 40?

    -   Note that many will fail/timeout or won't be useable
    -   Find links manually
    -   Create sam-dump commands from Rscript file
    -   Download to airlock drive

-   [ ] Use scripts from [sam2conseq](https://github.com/PoonLab/sam2conseq/blob/master/sam2conseq.py) to find consensus sequences for sam files

    -   Add these to the files to be sent into pangolin
    -   Potentially, find the reported sequence (as submitted to NCBI) instead

-   [ ] Ensure pangolin is working as expected

    -   Try it out on Rei?

-   [ ] Polish the visualizations of uncertainty

    -   Create new ones?

-   [ ] New folder: covid\_src

    -   Populate with existing source-able scripts, ensure outputs are consistent
    -   Add error-checking steps

-   [ ] Make all remaining scripts source-able, with appropriate outputs

-   [ ] Write pipeline script to source the files

    -   Does not need to include downloading scripts
    -   May be best to separate out parsing scripts as well (even with speedups, 30 sam files will take several days)
