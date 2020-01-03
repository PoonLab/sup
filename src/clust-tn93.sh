#!/bin/sh

x=$(ls seqs/seqs-prm-*.fasta)

for i in $x  
do
  tn93 -t 1 -f csv -o $i.out -q $i &
done


