#!/bin/sh

files=`ls data/sam/`

for file in $files; do
    src/parse-sam-c/sam2aln -f data/sam/$file -t 6
done

mkdir -p data/output
mv *.csv data/output
