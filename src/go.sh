#!/bin/sh

./clean.sh

# Generate the true phylogeny.
# Sequences saved in "seqs/sim.fasta"
# Tree saved in "trees/sim.nwk"
Rscript simul-seq.R

# Add uncertainty to the true sequences:
Rscript add-uncertainty.R

# Draw from the probabilistic sequences (Monte Carlo):
Rscript draw-tip-seqs.R

# Reconstruct phylogeny using FastTree:
./calc-tree.sh

# Compare the reconstruction
# to the true phylogeny:
Rscript ss.R

