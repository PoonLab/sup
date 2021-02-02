# /bin/sh, but not actually so no bang

# must be in the main fasta files directory

# TODO: Make this run the commands, not just create them
    # Check for existence first
    # associate the accession number with the run number
        # will be useful for posterity
# Give full path to binary to R, export doesn't work
Rscript ~OneDriveUWO/0postdoc/sup/create_sam_commands.R

# Check if the files have the fields I need,
# organize them into paired, unpaired, or badsam
Rscript sam_validator.R --tidyup

# Create the associated fasta file (requires samtools)
Rscript fasta_maker.R paired
Rscript fasta_maker.R unpaired

# Move files to Rei (too long on my computer)
scp -r paired dbecker7@rei:run_sam
scp -r unpaired dbecker7@rei:run_sam

# Parse the sam files
# Run on Rei
# TODO: take input/output directories as arguments
    # Stop outputting so much to nohup.out
    # Ensure parse.sam is updated on rei (maybe with git)
# nohup Rscript open_seSAMe_dir.R &

scp -r dbecker7@rei:run_sam/paired_DONE .
scp -r dbecker7@rei:run_sam/unpaired_DONE .




