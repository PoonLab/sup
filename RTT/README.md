# Root-To-Tip Regression for SARS-CoV-2 with Artificial Uncertainty

**Goal:** Calculate root-to-tip regression estimate for the TMRCA, then add artificial uncertainty at levels similar to our downloaded samples and observe the effect on the TMRCA.

Since TMRCA is just linear regression, there will be a standard error estimate for the intercept. Our analysis should show that this standard error is an underestimate of the true standard error. We also expect bias in the estimate itself. 

Research steps:

- [x] Collect the data (subsample according to world active case counts **or** uniformly across time for a manageable set)
	- [x] Ensure that the sampled data have associated SAM files and that these SAM files are valid.
- [x] Fit a phylogenetic tree to the raw data.
- [x] Find the ~~intercept~~ slope and its standard error from RTT regression
- [ ] Create N multi-fasta files from the uncertainty matrix tgat came from the SAM files
- [ ] For each multifasta, align to reference and fit a time tree using `treetime`
	- [ ] Record the slopes and their standard errors from each tree
- [ ] Summarise the difference between the slope and SE from conseqs versus the reslopes from the resamples (and possibly their SE's)
	- Possible summaries: violin plot of the slopes from resamples with a horizontal line for the original slope 
	- Single point with whiskers for standard errors for each resampled slope (with arbitrary index **or** order statistic on the x axis) and the slope and SE of the conseqs as a horizontal line and transparent box, respectively.

## Data

On June 23, 2021, we downloaded *all* data from:

https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&VirusLineage_ss=Severe%20acute%20respiratory%20syndrome%20coronavirus%202,%20taxid:2697049&HostLineage_ss=Homo%20(humans),%20taxid:9605&HostLineage_ss=Homo%20sapiens%20(human),%20taxid:9606&Completeness_s=complete

This data had the following information in the FASTA definition line:

    Accession | SRA Accession | Collection Date | Geo Location | Pangolin | Host | Isolation Source | Length | Assembly

These data were first filtered so that all entries had an associated SRA number (indicating that they were submitted to NCBI and thus the short run SAM files were available). We then sampled 5 genomes from each week and downloaded these sequences. We only used 2 per week, but 5 were downloaded to ensure the existence of the CIGAR string which is necessary for the construction of the uncertainty matrix but was missing from many of the files. 








