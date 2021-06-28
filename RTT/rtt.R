# Root to Tip Regression
library(ape)
library(dplyr)

# Load in tree ----------------------------------
tree <- read.tree("raw_seq.nwk")
tree$tip.label
# Tip labels only have accession for some reason
# Need to re-join with descriptions


# Load FASTA definition line --------------------
desc <- read.csv("sequences_descr_wk.csv")
desc$Accession <- trimws(desc$Accession)

# Sanity check
all(tree$tip.label %in% desc$Accession)

tips <- data.frame(label = tree$tip.label) %>%
	mutate(order = 1:n()) %>% 
	left_join(select(desc, Accession, date), by = c("label" = "Accession"))
tips

