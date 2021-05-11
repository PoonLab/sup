# Source Code for Simulated Uncertainty

To run the main analysis: `./run_multi.sh`

The general algorithm:

1. Simulate a root sequence
2. Simulates an evolutionary tree for that sequence
3. Generates uncertainty around the tree nodes
4. For each set of generated uncertain nodes, fit a tree
	- Also, calculates within-tree TN93
5. Summarizes the differences between inferred trees
	- Also difference between "certain" tree and inferred tree
6. Makes pretty plots

## Dependencies

- raxml (https://cme.h-its.org/exelixis/web/software/raxml)
- tn93 (https://github.com/veg/tn93)

## Files in this directory 

Files listed in order of evaluation.

- clean.sh: removes all tree files
- check-R-packages.R: checks for / installs R packages
	- I had a problem with Rscript using a different directory.
- simul-seq.R: Simulates the phylogeny
	- Relies on prm.csv for it's values
	- Outputs a series of trees
- RAXML is called from within run-multi
	- Outputs tree-raxml-*-mc.out (no time, unrooted)
- run-unit.sh
	- infr-tree-raxml.sh: inputs tree-raxml-*-mc.out
		- Outputs rooted time-trees
	- dist-tn93.sh: calculates TN93 distances for nodes WITHIN each tree
	- dist-calc.R: calculates TN93 distance statistics
		- formats them in a nice list
- analysis.R
















