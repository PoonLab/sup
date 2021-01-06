# sup - dev branch (short for devan, not development)

Sequencing Uncertainty Propagation

## Objectives

Propagate sequencing uncertainty in phylogenetic analysis, with an application to SARS-CoV-2 lineage assignment.

Sequencing is a multi-step process which is prone to errors. If the output for the sequencing of a biological sample gives `ATTGCTATGC`, what is the error probability associated with this result, at each position (for example, what is the probability that the second base is indeed a `T`)? How can we propagate this uncertainty in downstream phylogenic analysis? How variable are our estimates?

## Directory Structure

- `data`: Data related to sequencing uncertainty.
  - `seqs`: fasta sequence data for HIV data.
  - `zanini`: HIV patient data for Zanini et. al (2015). Population Genomics of Intrapatient HIV-1 Evolution. doi: https://doi.org/10.7554/eLife.11282.
  - `parsed_covid`: sequence uncertainty matrices (S) for SARS-CoV-2 data (raw sam files are too for GitHub, but filenames correspond to [Sequence Read Archive](https://www.ncbi.nlm.nih.gov/sra) runs).
  - `sampled_covid`: fasta files generated from `parsed_covid` 
- `doc`: Documentation of the methods.
- `src`: source code for uncertainty propagation analysis, based on simulated sequences with beta-distributed uncertainty at each site. 
  - `run-multi.sh` is where the magic happens.
  - `prm.csv` includes the parameters for the simulation and analysis.
- `reads-seq-err`: Estimation of the sequencing error of DNA fragment by Illumina instruments using simulations from from the software *InSilicoSeq*.
- `sung`: source code for `sung` (Sequence UNcertainty Generation) R package.
  - Actual R package is in `Rpkg` sub-directory.
- `covid`: source code for analysis of SARS-CoV-2 sam files. 
  - `open_seSAMe_dir.R` looks for `.sam` files in a directory and runs `parse-sam.r` (modified version of the one found in `sung`), then saves the results to `data/parsed_sam`. 
  - `sample_S.R` draws samples of sequences from the uncertainty matrices (S) in `parsed-covid`. These sequences can then be fed into [pangolin](https://github.com/cov-lineages/pangolin) using the codes found in `shell_commands.txt` 
    - TODO: streamline the analisys so that I don't need to switch between R and terminal)
  - `pangolineages.r` sends fasta from `sampled_sam` and 

(Note: I recognize that "covid" is the disease, "SARS-CoV-2" is the virus, but covid is easier to type.)

## Running scripts

### Main simulation analysis

`run-multi.sh` performs uncertainty analyses with different beta distribution shape values (i.e., uncertainty levels) as defined in `prm-btshp.csv`.

The type of sequence and phylogeny simulated is defined in `prm.csv`.

