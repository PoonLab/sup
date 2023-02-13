#!/bin/sh

make
files=`ls ./data`

for entry in $files; do
  ./sam2aln -f data/$entry -t 3
done
