#!/bin/zsh

MSNAME=sup-ms

pdflatex $MSNAME.tex
bibtex $MSNAME
pdflatex $MSNAME.tex
pdflatex $MSNAME.tex
	
