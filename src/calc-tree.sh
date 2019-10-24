#!/bin/sh

# Reconstruct the phylogeny from 
# sampled probabilistic sequences

# How many Monte Carlo samples:
N=$(ls -l seqs/seq*fasta | wc -l)
echo Calculating phylogeny for $N trees...

# For each set of sampled tips, 
# reconstruct the phylogeny:
for i in $(seq 1 $N) 
do
 FastTree -quiet -nt seqs/seqs-mc-$i.fasta > trees/tree-mc-$i.out
done

echo Reconstruction of all MC trees with FastTree done.

