# Workflow Outline

This is generally how I'm approaching the steps. I'm writing this for my own benefit, since my workflow recently changed and I need to remind myself what I'm doing.

## Step 1: Access the Files

### SAM Files

1. On [NCBI Short Run Archive](https://www.ncbi.nlm.nih.gov/sra/?term=txid2697049%5BOrganism:noexp%5D%20NOT%200[Mbases), select the checkboxes for an arbitrarily chosen (by me) selection of files.
    - Search parameters: Illumina (best chance of having Cigar strings), genome, bam, optionally paired/unpaired.
2. Select "Send to:" at the top of the page, choose File, choose RunInfo.
    - Downloads a csv with all of the information about the selected runs.
    - I save this file to my external hard drive.
3. Prepare folder structure.
    - Folders for paired and unpaired
4. Run Downloader.R in the same working directory as the RunInfo csv and paired and unpaired folders.
    - This script will loop through the rows, attempt to download the csv, and report any errors that arise (written in a copy of the RunInfo csv).
    - This script has completely filled my external drive.

The files will be neatly organized, with any errors deleted.

I went a step further and put the fasta files into one folder, and split the sam files over multiple folders that are about 10GB in size. This makes it easier to send them to Rei, and on Rei the folders can be parsed one at a time.

### Called Sequences (WIP)

Using samconseq?
Modify sample_S to use the max probability?



