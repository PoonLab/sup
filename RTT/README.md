# Root-To-Tip Regression for SARS-CoV-2 with Artificial Uncertainty

**Goal:** Calculate root-to-tip regression estimate for the TMRCA, then add artificial uncertainty at levels similar to our downloaded samples and observe the effect on the TMRCA.

Since TMRCA is just linear regression, there will be a standard error estimate for the intercept. Our analysis should show that this standard error is an underestimate of the true standard error. We also expect bias in the estimate itself. 

Research steps:

- Collect the data (subsample according to world active case counts for a manageable set)
- Fit a phylogenetic tree to the raw data.
- Find the intercept from RTT regression
- For a range of uncertainty values, do:
    - For N times, do:
        - Add uncertainty to all of the sequences according to the beta distribution
        - Re-sample from the uncertainty 
        - Fit a phylogenetic tree to the new samples
        - Calculate the intercept (and it's standard error) in the root-to-tip regression



## Data

On June 23, 2021, we downloaded *all* data from:

https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&VirusLineage_ss=Severe%20acute%20respiratory%20syndrome%20coronavirus%202,%20taxid:2697049&HostLineage_ss=Homo%20(humans),%20taxid:9605&HostLineage_ss=Homo%20sapiens%20(human),%20taxid:9606&Completeness_s=complete

This data had the following information in the FASTA definition line:

    Accession | SRA Accession | Collection Date | Geo Location | Pangolin | Host | Isolation Source | Length | Assembly

The data will be subsampled so that the phylogenetic tree can be constructed in a reasonable amount of time (there are somewhere between 5 to 20 parameter sets and N replications for each parameter set).

For my future reference, the description file was created with the following shell commands:


## Adding Artificial Uncertainty

Using the uncertainty sequences that we have already calculated, we can estimate the average per-base error rate for several files. The parameters of the beta distribution will be based on these error rates. The process for this is outlined in sup/ms. 







