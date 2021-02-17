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

`data`: Data related to sequencing uncertainty

`doc`: Documentation of the methods

`src`: source code for uncertainty propagation analysis

`reads-seq-err`: Estimation of the sequencing error of DNA fragment by Illumina instruments using simulations from from the software *InSilicoSeq*.


## Running scripts

### Main analysis

`run-multi.sh` performs uncertainty analyses with different beta distribution shape values (i.e., uncertainty levels) as defined in `prm-btshp.csv`.

The type of sequence and phylogeny simulated is defined in `prm.csv`.


### Fragment sequencing eror
