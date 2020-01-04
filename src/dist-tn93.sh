#!/bin/sh

# CALCULATE THE TN93 DISTANCE
# FOR ALL SEQUENCES WITHIN 
# EACH TREE (MCxPRMSET)

x=$(ls seqs/seqs-prm-*.fasta)

for i in $x  
do
  tn93 -t 1 -f csv -o $i.out -q $i &
done


