#!/bin/sh

files=`ls data/sam/`

for file in $files; do
    src/parse-sam-c/sam2aln data/sam/$file 6
done

mkdir -p data/output
mv *.csv data/output
