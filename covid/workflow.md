# Workflow Outline

This is generally how I'm approaching the steps. I'm writing this for my own benefit, since my workflow recently changed and I need to remind myself what I'm doing.


## Step 1: Access the SAM Files

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


# Step 2: Parse the files

1. Send the files to Rei
    - Especially in a chunk. The folders labelled "blockp" and "blocku" are folders containing 10 gigabytes of paired and unpaired sam files, respectively.
2. Run parse.sam, using a script that searches for completed files then loops through incomplete ones.
    - nohup Rscript open_seSAMe_dir.R blocku > nohup1.out &
    - Files are placed in samDONE, regardless of where they came from. Repeats are replaced.
3. The resulting RDS files are small enough to be stored on github.


# Step 3: Sample from the files

1. sample_S_dir.R looks through samDONE and samples from anything without a sample.
    - The final lines of the script generate the pangolin calls
    - This script also calculates the consensus sequence (putting Ns if there are fewer than 10 observations at a given location)

TODO: clean up the script, separate the generation of pangolin calls
    

# Step 4: pangolin calls

1. conda activate pangolin; copy and paste the pangolin calls from sample_S_dir.R


# Step 5: pangolin_results_dir.R (WIP)

There are multiple attempts at pretty plots showing the lack of concordance.

TODO: Choose/make better ones.


# Step 6: There is no step 6

Drink a beer, maybe?




