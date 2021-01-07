# sup - dev branch (dev is short for devan, not development)

Sequencing Uncertainty Propagation

## Objectives/Overview

Propagate sequencing uncertainty in phylogenetic analysis, with an application to SARS-CoV-2 lineage assignment.

Sequencing is a multi-step process which is prone to errors. If the output for the sequencing of a biological sample gives `ATTGCTATGC`, what is the error probability associated with this result, at each position (for example, what is the probability that the second base is indeed a `T`)? How can we propagate this uncertainty in downstream phylogenic analysis?

Sequence uncertainty can be obtained either from SAM files or from FASTQ:

1. Raw short read files (e.g. SAM)
  - Sequences are read little bits at a time. Each read is recorded, and stored in a file (along with it's alignment to a reference sequence).
  - The sequence reads include information on the read quality. This is encoded as a Phred score, which uses unicode characters to encode uncertainty about a given base call.
  - The called base is assigned a probability according to the Phred score. The remaining probability is assigned equally to the remaining possible bases.
    - For example, if P(T) = 0.7 from the Phred score, then P(A) = P(C) = P(G) = 0.1.
  - Given many short reads, the probabilities for each base pair at each site are added.
    - If, say, T is always 100% certain, then the the sum of the probabilities will represent the number of reads.
2. FASTQ Files
  - SAM files aren't always practical (or available), so often FASTQ files are used instead.
  - These contain a Phred score for each base call. This phred score is calculated similarly to the process for SAM files, but only the uncertainty for the most probable base is reported.
3. FASTA files
  - Contain no information about sequence uncertainty.

In this study, we show that using sequences from FASTA files leads to underestimation of the variance, which has ripple effects throughout the rest of the analysis. We demonstrate techniques for propagating sequence uncertainty into further analysis, but this inevitably comes at the expense of computation time.

We include an application to SARS-CoV-2 data. Using a collection of SAM files for the SARS-CoV-2 virus from [NCBI's short read archive](https://www.ncbi.nlm.nih.gov/sra), we produce uncertainty estimates at each site on the genome. Instead of choosing the most likely base at each site to create a single nucleotide (as in FASTA files), we sample a collection of sequences based on the uncertainty. [Pangolin](https://github.com/cov-lineages/pangolin) is used to assign each sequence to a lineage.

## Directory/Analysis Structure

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
  - Written by David Champredon and Art Poon.
- `covid`: source code for cluster allocation analysis of SARS-CoV-2 sam files. (Note: I recognize that "covid" is the disease and "SARS-CoV-2" is the virus, but covid is easier to type.)
  - `open_seSAMe_dir.R` looks for `.sam` files in a directory and runs `parse-sam.r` (modified version of the one found in `sung`), then saves the results to `data/parsed_sam`.
    - Note: `parse-sam.r` takes about 5 hours on Rei for a 200MB SAM file.
  - `sample_S.R` draws samples of sequences from the uncertainty matrices (S) in `parsed-covid`. These sequences can then be fed into [pangolin](https://github.com/cov-lineages/pangolin) using the codes found in `shell_commands.txt`
    - TODO: streamline the analisys so that I don't need to switch between R and terminal.
    - TODO: samples from S should be based on posterior distribution (Neg Binom) to incorporate number of reads at each site.
  - `pangolin_results.r` analyses output from pangolin, especially with respect to variance.
- `gsun`: a few files for uncertainty generation testing
- `misc`: a few files for calculating testing RF distance
- `ms`: start of the paper






