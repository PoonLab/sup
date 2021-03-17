### Specific to this branch - con branch (Following naming convention)
Primarily working on rewriting the parse-sam.R code to increase speed and memory management. This makes 4 general changes... 
1. Initially read in files by chunks, using `stringr::read_lines_chunked()` (see `parse.sam()`)
2. Use data tables to store results (including data read from sam) and use data table row operations whenever possible.
3. Vecotrize functions where possible (see - `parse.sam.line()`, `apply.cigar() )
4. Make use of parallel capabilities through `parallel::mclapply()` for file reading and vectorized functions.

# Original README
# sup
Sequencing Uncertainty Propagation

## Objectives

Propagate sequencing uncertainty in phylogenetic analysis.

Sequencing is a multi-step process which is prone to errors. If the output for the sequencing of a biological sample gives `ATTGCTATGC`, what is the error probability associated with this result, at each position (for example, what is the probability that the second base is indeed a `T`)? How can we propagate this uncertainty in downstream phylogenic analysis?

## Directories


1. Raw short read files (e.g. SAM)
    - Sequences are read little bits at a time. Each read is recorded and stored in a file (along with it's alignment to a reference sequence).
    - The sequence reads include information on the read quality. This is encoded as a Phred score, which uses unicode characters to encode uncertainty about a given base call.
    - The called base is assigned a probability according to the Phred score. The remaining probability is assigned equally to the remaining possible bases.
        - For example, if P(T) = 0.7 from the Phred score, then P(A) = P(C) = P(G) = 0.1.
    - Given many short reads, the probabilities for each base pair at each site are added together.
        - If, say, T is always 100% certain, then the the sum of the probabilities will represent the number of reads.
2. FASTQ Files
    - SAM files aren't always practical (or available), so often FASTQ files are used instead.
    - These contain a Phred score for each base call. This phred score is calculated similarly to the process for SAM files, but only the uncertainty for the most probable base is reported.
3. FASTA files
    - Contain no information about sequence uncertainty.


`doc`: Documentation of the methods

`src`: source code for uncertainty propagation analysis

`reads-seq-err`: Estimation of the sequencing error of DNA fragment by Illumina instruments using simulations from from the software *InSilicoSeq*.

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

# Roadmap to the Paper

- WRITE
    - Currently making LOs and first drafts by topics, will put it all together soon.
    - Intro
        - [ ] Intro paragraphss
        - [ ] Relevant literature
        - [ ] Sources of uncertainty
        - [ ] Description of SAM files (with and without paired reads)
        - [ ] Description of InSilicoSeq (and other things from David's work)
    - Methods
        - [ ] Uncertainty Matrices (and sequence-level uncertainty, similar to likelihoods; largely re-used from David's work)
        - [ ] Sampling from normalized probabilities or from multinomial posteriors (justification)
        - [ ] Sequence-level uncertainty measures, and their use in analysis of genetic data
    - Application
        - [ ] Description of Covid Data
        - [ ] Example of an uncertainty matrix from actual covid data
        - [ ] Pangolin description (incl. bootstrap support)
        - [ ] Pangolin results summary/vis
    - Conclusions
        - [ ] For phylogenies in general (Lean heavily on David's work)
        - [ ] For Covid data
- `parse.sam()` is slow and unreasonable for our purposes. 
    - Me, Connor, and Gopi are all working on methods to make it faster. 
- [x] Visualizations of uncertainty matrices (`covid/unc_vis.R`)
    - Mostly done, but needs some annotations.
- [ ] Fixing the sampling code (`covid/sample_S`)
    - Changed to multinomial posterior and now resultant samples are all 100% certain according to pangolin. 
        - Possible fixes: Higher coverage cutoff (so that there are more Ns)
    - [ ] Compare samples with/without posterior sampling
- [ ] Visualizations of Pangolin calls (`covid/pangolin_results_dir.R`)
    - Needs a lot of work.
        - Want: Find genetic distance between conseq and other sequences (to validate my sampling)






