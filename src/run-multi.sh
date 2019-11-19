#!/bin/sh

./clean.sh


# Generate the true phylogeny.
# Sequences saved in "seqs/sim.fasta"
# Tree saved in "trees/sim.nwk"
Rscript simul-seq.R

for i in {1..3}
do
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
  Rscript ss.R $i
done




