#!/bin/sh

make
files=`ls ./data`

for entry in $files; do
  ./sam2main data/$entry 3
done
