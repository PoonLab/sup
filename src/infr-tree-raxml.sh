#!/bin/sh

# Reconstruct the phylogeny from 
# sampled probabilistic sequences

# RAxML documentation:
# https://cme.h-its.org/exelixis/web/software/raxml/index.html

figgy=`whereis figlet`
figgylen=${#figgy}
if [ $figgylen > 10 ]; then 
    figlet Analysis.R
fi

PRMSET=$1

# How many Monte Carlo samples:
N=$(ls -l seqs/seqs-prm-$PRMSET*.fasta | wc -l)
echo "Calculating phylogeny for $N trees (raxml)..."

# For each set of sampled tips, 
# reconstruct the phylogeny:
RAXML_OPT="-m GTRGAMMA -p 12345 --JC69 --silent"

# Trees infered from drawing
# probabilistic sequences:
TREE_NAME=tree-raxml-prm-$PRMSET-mc

for i in $(seq 1 $N) 
do
  # Phylogeny inference:
  raxmlHPC $RAXML_OPT -s seqs/seqs-prm-$PRMSET-mc-$i.fasta -n $TREE_NAME-$i.out > trees/RAxML-prm-$PRMSET-mc-$i.out
done

echo "Unrooting ..."
# Remove the `:0.0` string that 
# symbolizes rooting in RAxML tree:
for i in $(seq 1 $N) 
do
  sed -i -e 's/:0.0;/;/g' RAxML_bestTree.$TREE_NAME-$i.out
done

# RAxML can only output in the local directory,
# so move output files manually:
mv RAxML_*$TREE_NAME*.out* trees/
rm trees/*.out-e #TODO: fix that

echo "Reconstruction of all MC trees with RAxML done."
echo "RAxML outputs saved in trees/tree-mc-raxmli.out."

