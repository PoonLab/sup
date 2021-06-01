#!/bin/sh

make
files=`ls ./data`

for entry in $files; do
  ./sam2aln data/$entry 3
done
