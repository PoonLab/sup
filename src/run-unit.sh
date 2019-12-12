#!/bin/sh

i=$1

echo "Processing parameter set #$i ..."

# Add uncertainty to the true sequences:
Rscript add-uncertainty.R $i

# Draw from the probabilistic sequences (Monte Carlo):
Rscript draw-tip-seqs.R $i

# Reconstruct phylogeny using FastTree:
#./infr-tree-fasttree.sh > fasttree.out

# Reconstruct phylogeny using RAxML:
./infr-tree-raxml.sh $i > raxml.out

# Compare the reconstruction
# to the true phylogeny:
Rscript dist-calc.R $i
