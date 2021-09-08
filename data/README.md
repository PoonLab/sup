# Data folder

## Sub-folders

- unc_covid:
	- The RDS files that result from running parse-sam-c on the SAM files
		- Also results in "*_insertions.RDS" when there are insertions present.
	- NOT tracked by git, since these files are large, they change, and there are a lot of them
- sampled_covid:
	- The result of running `bash covid/src_sample-S.R`, generally with `-N 1000`
		- `-d` in the call above results in dirichlet sampling, which also appends a `_d` to the filenames.
	- Results in the csv files in pangolineages.
- pangolineages:
	- output of pangolin, labelled with the accession number and a d for dirichlet sampling
- pangordlineages:
	- output of pangolin, but for ordered likelihood analysis
- zanini:
	- example data from Zanini et. al. 2015
	- R script (TODO: explain this code)
