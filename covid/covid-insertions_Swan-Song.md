# Swan Song for Insertions

A great deal of thought and effort has been applied to the problem of incorporating insertions into the re-sampling algorithm. Art Poon, Gopi Gugan, and Devan Becker (me) have all made excellent efforts in this direction.

Alas, the problem still requires a much more sophisticated solution to be viable, and thus is being removed from the current codebase (and manuscript) and will be added to the future work section.

# Reason

The handling of the insertions was based on our naive - but currently necessary - assumption of independence between sites. In other words, the insertions are being added without regard for what's happening at other sites.

We note that all of the insertions are biologically plausible. If we sample an insertion for a given site, that insertion is one that was actually recorded in at least one of the short reads. Insertions of length, say, 2, were always treated as a length-2 insertion; we never just used one of the two inserted base pairs.

Despite our best efforts, the re-sampled genomes with the insertions resulted in genomes that the Pangolin system was unable to map. This is likely due to the fact that insertions in the SARS-CoV-2 genome are rare events and were always preceeded by a specific SNP or set of SNPs elsewhere in the genome. When Pangolin tries to align genomes with insertions to the other observed genomes, it is missing these other SNPs and thus cannot map the genome anywhere in the tree. 

We note that the tree used by Pangolin is based on the assumption that the consensus sequences are known deterministically, which means that there are almost certainly sequences that have errors in their tree. It is possible that more of our resampled sequences would be mapped to the tree if these errors did not exist. The purpose of this work is to expose these errors, so it is perhaps ironic that our need to expose a phenomenon is hindered by that phenomenon not having been exposed yet.

# Future Work

To rectify these issues, we would need to develop an algorithm that takes into account correlation between position. Specifically, any re-sampling of the genome would need to condition on all other locations. At the bare minimum, we believe that any resampling must allow for mutations that always occur with another mutation (although this relationship need not be commutative - an insertion may only occur with a given SNp, but that SNP may occur without the insertion). This is feasible, but would require careful consideration of the algorithm and the required data structure.

An alternative approach may eschew re-sampling altogether. For instance, a list of all biologically feasible genomes could be constructed based on all possible combinations of the short reads in a SAM file, and then these could be paired with their genome likelihoods. Determining which of the possible combinations are biologically feasible would be a non-trivial task, and there are likely millions - if not billions or even trillions - of biologically possible combinations of the short reads in a SAM file. 

It is also likely that our insertions were not being added at the correct rate. In order to sample the presence of an insertion, we needed to know both the numerator (the number of observed insertions at that location) and the denominator (the number of reads at the position *before* the insertion). Our algorithm for parsing the SAM files was not prepared to find the denominator, and thus this was estimated based on the coverage (by construction, the sum of each column of the uncertainty matrix represents the number of reads at that position). Our approach treats paired reads as two halves of a single observation, but the number of insertions does not take paired reads into account. This disconnect means that the numerator is a count of insertions (where pairs a full observation) but the denominator is a count of observations (where pairs are half of an observation). Furthermore, the coverage at the position before the insertion is not a perfect measure of the number of potential reads at the insertion position. Any continuation of this work would need to address these issues.
