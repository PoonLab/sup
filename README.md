# sup - Sequencing Uncertainty Propagation

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
    - However, uncertainty can be generated around FASTA to investigate how sensitive the methods are to the assumption of certainty in the conseq.


- `drafts`: Skelton.Rmd is the draft for the final paper.
- `covid`: source code for cluster allocation analysis of SARS-CoV-2 sam files. (Note: I recognize that "covid" is the disease and "SARS-CoV-2" is the virus, but covid is easier to type.)
    - `open_seSAMe_dir.R` looks for `.sam` files in a directory and runs `parse-sam.r` (modified version of the one found in `sung`), then saves the results to `data/parsed_sam`.
        - Note: `parse-sam.r` takes about 5 hours on Rei for a 200MB SAM file.
        - **Deprecated**, use `parse-sam-c` in src instead.
    - `sample_S.R` draws samples of sequences from the uncertainty matrices (S) in `parsed-covid`. To be fed into [pangolin](https://github.com/cov-lineages/pangolin).
    - `ordered-lik.R` investigates calculating the *n* most likely sequences to see at which likelihood value the results of the analysis changes.
    - `pangolin_results.r` analyses output from pangolin, especially with respect to variance.
- `src`: source code for uncertainty propagation phylogeny analysis, based on simulated sequences with beta-distributed uncertainty at each site.
    - `run-multi.sh` is where the magic happens.
	- Generates a phylogeny, generates uncertainty on top of this phylogeny, estimates the trees, then calculates the distance between each tree and the tree without uncertainty.
    - `prm.csv` includes the parameters for the simulation and analysis.
- `data`: Data related to sequencing uncertainty.
    - `seqs`: fasta sequence data for HIV data.
    - `zanini`: HIV patient data for Zanini et. al (2015). Population Genomics of Intrapatient HIV-1 Evolution. doi: https://doi.org/10.7554/eLife.11282.
    	- `multinomial-zanini.R` is the workhorse
    - `parsed_covid`: sequence uncertainty matrices (S) for SARS-CoV-2 data (raw sam files are too for GitHub, but filenames correspond to [Sequence Read Archive](https://www.ncbi.nlm.nih.gov/sra) runs).
    - `sampled_covid`: fasta files generated from `parsed_covid`
    - `ord_covid`: 
- `reads-seq-err`: Estimation of the sequencing error of DNA fragment by Illumina instruments using simulations from from the software *InSilicoSeq*.
- `doc`: Documentation of the methods.
- `sung`: source code for `sung` (Sequence UNcertainty Generation) R package.
    - Written by David Champredon and Art Poon.
- `gsun`: a few files for uncertainty generation testing
	- Early draft files
- `misc`: a few files for calculating testing RF distance
- `ms`: start of the paper

# Roadmap to the Paper

- WRITE
    - Intro
        - [x] Intro paragraphss
        - [x] Relevant literature
        - [x] Sources of uncertainty
        - [x] Description of SAM files (with and without paired reads)
        - [x] Description of InSilicoSeq (and other things from David's work)
        - [ ] More citations to support opinions
        - [ ] Revise & rewrite
    - Methods
        - [x] Uncertainty Matrices (and sequence-level uncertainty, similar to likelihoods; largely re-used from David's work)
        - [x] Sampling from normalized probabilities or from multinomial posteriors (justification)
        - [x] Sequence-level uncertainty measures, and their use in analysis of genetic data
        - [ ] Beta distribution to generate uncertainty
    - Application
        - [x] Description of Covid Data
        - [x] Pangolin description (incl. support)
        - [x] Pangolin results summary/vis
        - [x] Intro to David's work
        - [x] Goal of David's work
        - [ ] Results from David's work (after re-running)
        - [ ] Takeaways from David's work
    - Conclusions
        - [ ] For phylogenies in general (lean heavily on David's work)
        - [x] For Covid data
		- [x] Other uses: instead of a canonical representation, the probability sequence could be calculated based on the estimated population frequency of alleles at polymorphism sites.

### Nice-to-Haves

- Comparison of sampling versus ordered likelihoods
- Quantification of the actual levels of uncertainty (e.g. histogram).
- 




