# Handling insertions in the construction of the uncertainty matrix

Let's begin with the following toy example. To make things explicit, I am first going to represent it as an aligned SAM file where all of the short reads are aligned to each other, which results in padding of the reference sequence (in practice, short reads are only aligned to the reference; multiple alignment would need to be done with respect to hundreds of thousands -or millions - of short reads and thus I do not consider that). 

The stars in the position number (column titles) represent padding in the reference due to alignment of insertions in reads. Blanks in the columns mean that position was not read, whereas gaps (`-`) represent deletions. If there are two characters in the entry, the second digit represents the phred score (`%` means 50%, `)` represents 20% probability of error). For simplicity, if there is no second digit then there is essentially no chance of error (Phred score of `I`). 

| Read | 1 | 2 | 3 | 4  | *  | *  | 5 | 6 | 7 |
|------|---|---|---|----|----|----|---|---|---|
| Ref  | A | T | C | G  | -  | -  | C | T | A |
| 1    |   |   | C | -  | A) | -  | C | T |   |
| 2    | A | T | C | G) |    |    |   |   |   |
| 3    |   |   |   |    |    |    | C | T | G |
| 4    |   |   | C | G  | T% | T% | C |   |   |
| 4    |   |   | C | G  | A% | T% | C |   |   |
| 5    |   | T | - | G  | A  | -  | C | T |   |
| 6    |   |   |   | G% | -  | T  | C | - | A |
| 7    |   |   |   |    | A  | -  | C | - | A |
| 8    |   |   | C | G  | -  | -  | C | T | A |

The two padded positions (pads) are designed to represent a "true" insertion of AT, but some reads saw a deletion of the A and some saw a deletion of the T. Some other notes:

- Read 2 has coverage up to position 4, but not coverage of the pads.
- Read 3 has coverage after the pads.
- Read 7 has coverage of the pads, but not at position 4.

In an ideal world, the table above would represent our data. However, the table below is what we actually get:

| READ | POS | CIGAR      | SEQ    | Q     |
|------|-----|------------|--------|-------|
| 1    | 3   | 1M1D1I2M   | CACT   | I)II  |
| 2    | 1   | 4M         | ATCG   | III)  |
| 3    | 5   | 3M         | CTG    | III   |
| 4    | 3   | 2M2I1M     | CGTTC  | II%%I |
| 4    | 3   | 2M2I1M     | CGATC  | II%%I |
| 5    | 2   | 1M1D1M1I2M | TGACT  | IIIII |
| 6    | 4   | 1M1I1M1D1M | GTCA   | %III  |
| 7    | 4   | 1I1M1D1M   | ACA    | III   |
| 8    | 3   | 5M         | CCGTA  | IIIII |

Note that these two tables are *not* necessarily interchangable. In order to re-construct the pads, the reads would need to be aligned to each other (the CIGAR is the alignment with the reference; note how read 8 appears to simply match the reference with no mention of insertions or deletions at position 4). This process is imperfect, and we may not always recover the correct insertion location relative to the other reads when there are multiple insertions in a row.

This table is used to construct the uncertainty matrix, which currently ignores the insertions and adds deletions as if they had an error probability of 0 (presumably, any detected base would have been reported and given a suitably large error probability and therefore we must be confident that there was no nucleotide to read). Because of the multiple alignment issue, the locations, error probabilities, and base calls of the insertions are recorded in a separate file rather than being incorporated into the uncertainty matrix. 

# Multi-stage sampling of insertions

To incorporate insertions of varying lengths, we first sample the number of insertions, then sample which nucleotides to insert based on insertion-specific uncertainty matrices. 

**Step 0:** Choose the denominator.

We want the denominator in the following steps to be indicative of the coverage at the insertion sites. However, we have no easy way of knowing how many reads would have had insertions at that position. Instead, it is estimated based on the coverage at the position before the insertions. 

We choose the position before somewhat arbitrarily. The example above shows reads that have coverage at position 4 but not in the pads, coverage at the pads but not position 4, and coverage at position 5 but not the pads. 

For the example above, the coverage at position 4 can be found as the sum of the phred scores, which is 4.2 (recalling that paired reads count as a full observation together; i.e., a half an observation apart). 

**Step 1:** Sample the number of insertions. 

The probability of inserting a single nucleotide can be be estimated as the number of single insertion events divided by the coverage at the position before the insertion. 

There are 4 single nucleotide insertions and the sum of their phred scores is 2.8, so the (un-normalized) rate of single insertions is 2.8/3.2 = 7/8. 

For two insertions, we cannot simply add the phred scores, so instead we first multiply them. For the first row of read 4 our example, this becomes 0.5 times 0.5, then multiplied by 0.5 again since this a paired read. The second read has the same phred scores, so the numerator here is 0.5^3 + 0.5^3 = 0.25.

For 0 insertions, there is a single read. However, there are no phred scores associated with this. To be consistent with our treatment of deletions, we assume that the error probabilities are 0 (which is also consistent with ignoring the deletions for the single nucleotide case). This leaves us with a numerator of 1.

To be consistent with the Dirichlet-multinomial sampling scheme, we can sample probabilities for each insertions according to a Dirichlet distribution with parameters 1/4.2, 0.25/4.2, and (7/8)/4.2. From these probabilities, we can sample the number of insertions. 

**Step 2:** Sample nucleotides from insertion-specific uncerainty matrices.

With the number of insertions determined, the insertion-specific uncertainty matrices can be created based on the insertions matching the criteria. For single insertions, the uncertainty matrix is as follows:

|   | 1   |
|---|-----|
| A | 2.8 |
| T | 1   |
| C | 0   |
| G | 0   |

Note that this does not (and cannot reasonably) account for the alignment of the reads in the two insertions.

For two insertions, the uncertainty matrix is:

|   | 1    | 2   |
|---|------|-----|
| A | 0.25 | 0   | 
| T | 0.25 | 0.5 | 
| C | 0    | 0   | 
| G | 0    | 0   | 

This matrix is based entirely on the paired read (Read 4), and so the columns would add to 1 if there were no error probabilty. However, both pairs of this read also had error probabilities of 0.5.

Note that the second column corresponds to always sampling a T. In general (and especially for positions with few insertion events), this process is amounts to inserting the observed insertions according to how often there were insertions.

Both of these matrices would be sampled from in the usual Dirichlet-multinomial method.

# Conclusion

This method is a reasonable approach given the data we have access to. We incorporate the uncertainty as best we can and the Dirichlet-multinomial distribution accounts for uncertainty due to coverage. 

The concessions we must make to account for the alignment to the reference - rather than to other short reads - highlights the importance of choosing a good reference sequence. If only there were a good clustering method that incorporated spatial information so that researchers could choose a reference sequence from the most likely cluster...


# Alternatives

- The numerators are based on phred scores, but could also be based on raw counts. This would lead to many improper fractions, but we have no requirement that the results add to 1 since they become parameters in the Dirichlet distribution.






