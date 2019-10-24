#!/bin/sh

N=$(ls -l seqs/seq*fasta | wc -l)
echo Calculating phylogeny for $N trees...

for i in $(seq 1 $N) 
do
 FastTree -quiet -nt seqs/seqs-$i.fasta > trees/tree-$i.out
done

echo FastTree done.

