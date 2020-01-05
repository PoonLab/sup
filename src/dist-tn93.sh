#!/bin/sh

# CALCULATE THE TN93 DISTANCE
# FOR ALL SEQUENCES WITHIN 
# EACH TREE (MCxPRMSET)

# To install tn93, see: https://github.com/veg/tn93

echo "Start TN93 distances calculations..."

x=$(ls seqs/seqs-prm-*.fasta)

for i in $x  
do
  tn93 -t 1 -f csv -o $i.out -q $i &
done

echo "TN93 distances calculated."

