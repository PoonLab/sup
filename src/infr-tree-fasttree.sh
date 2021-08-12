#!/bin/sh

# Reconstruct the phylogeny from 
# sampled probabilistic sequences

# FastTree documentation:
# http://www.microbesonline.org/fasttree/

# How many Monte Carlo samples:
N=$(ls -l seqs/seqs-mc*.fasta | wc -l)
echo Calculating phylogeny for $N trees fasttree...

# For each set of sampled tips, 
# reconstruct the phylogeny:

FT_OPT="-quiet -nt -spr 4 -mlacc 2 -slownni"

for i in $(seq 1 $N) 
do
 FastTree $FT_OPT seqs/seqs-mc-$i.fasta > trees/tree-mc-$i.nwk
done

echo "Reconstruction of all MC trees with FastTree done."
echo "FastTree outputs saved in trees/tree-mc-i.out."

