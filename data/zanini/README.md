# Data based on Zanini 2015

Main analysis: `invariant-prop.R`



## Files (in order of execution)

- utils.zanini.R: 
	- read_data: reads in all files of the form "act_p[patient]/[timepoint]" and "patient[patient]_day_[timepoint]", formats them nicely.
	- calc_entropy_one: sum(p * log_2(p))
	- calc_entropy: loops over positions and bases
	- digest: transfoms dataframes (?)
	- proba_by_base: uses digest, further transforms data and plots the error probability at each nucleotide location
	- plots: produces a bunch of plots
- multinomial-zanini.R 
	- calculates the entropy based on the multinomial distributionb
- invariant-prop: estimate the proportion of nucleotide locations which are highly unlikely to be errors






