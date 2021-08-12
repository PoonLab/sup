# SARS-CoV-2 Uncertainty Analysis

Calculate the sequence uncertainty matrix, sample from it (or calculate ordered likelihood), then feed it into pangolin.

- Main analysis (sampling): `src_sample-S.sh`
- Secondary analysis (ordered): `src_ordered-lik.sh`

## Dependencies

- pangoLEARN

## Files

Accessing SAM files:

- downloader.R !!!Not on GitHub!!!
	- After accessing NCBI SRA and recording accession, downloads files into a directory and records auxilliary information
- shell_commands.txt: a reminder of the commands I use (I'm new to this, okay?)
- samDONE: a list of SAM files that I'm waiting to parse, including ones with errors

Calculating the uncertainty matrix:

- open_seSAMe_dir.R
	- runs parse.sam on all files in a given directory
	- (run from command line)
	- open_seSAMe_single.R is the original file for testing
- parse-sam.r: contains functions (namely parse.sam) to covert SAM to S matrix

For Sampling analysis:

- src_sample-S.sh calls:
	- sample-S.R: samples from the uncertainty
		- Calls it once with overwrite = TRUE, then four more times
		- This uses less memory
	- pangolin for the 5000 sequences
	- Creates an Rmd with summary stats and plots

For Ordered likelihood analysis:

- src_ordered-lik.sh calls:
	- ordered-lik.R calculates the top $N most likely sequences, based on S
	- sends them through pangolin
	- Creates an Rmd with summary stats and plots 







