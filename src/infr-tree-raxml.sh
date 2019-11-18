#!/bin/sh

# Reconstruct the phylogeny from 
# sampled probabilistic sequences

# RAxML documentation:
# http://www.microbesonline.org/fasttree/

# How many Monte Carlo samples:
N=$(ls -l seqs/seqs-mc*.fasta | wc -l)
echo Calculating phylogeny for $N trees...

# For each set of sampled tips, 
# reconstruct the phylogeny:

RAXML_OPT="-m GTRGAMMA -p 12345 --JC69"

TREE_NAME="tree-raxml-mc"

for i in $(seq 1 $N) 
do
 raxmlHPC $RAXML_OPT -s seqs/seqs-mc-$i.fasta -n $TREE_NAME-$i.out
done

mv RAxML_*$TREE_NAME*.out trees/

echo "Reconstruction of all MC trees with RAxML done."
echo "RAxML outputs saved in trees/tree-mc-raxmli.out."

