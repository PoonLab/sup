#!/bin/sh

./clean.sh


# Generate the true phylogeny.
# Sequences saved in "seqs/sim.fasta"
# Tree saved in "trees/sim.nwk"
Rscript simul-seq.R

N=$(wc -l < prm-btshp.csv)

for i in $(seq 1 $N)  #TO DO: dont hardcode!
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

# Analyze the distances:
Rscript analysis.R

echo "  run-multi.sh completed. "

