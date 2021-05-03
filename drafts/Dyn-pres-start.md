# Propagating Sequence Uncertainty to SARS-CoV-2 Lineage Classification

Principles/Plan
- 4 page slideshow, more information than usual per slide (they can pause if they like)
- Motivation, Pangolin, Results


## Motivation

SRA: Read small parts of the sequence, higher variation but more reads

SAM file ---> Lose Uncertainty ---> Conseq ---> Pangolin ---> Lineage (D.7.2)

## Process

SAM ---> Sample 5000 ----> Run 5000 Pangolineage assignments
        (or top 5000)

## Results

- Lots of variation!!!
    - 100% boot probability
- Conseq almost never appears
- B.1.1.7 is always 100% (yay)




# Slide 1

SAM Files to Conseq to Analysis

- Some notes about the uncertainty



# Slide 2

SAM Files to Matrix to Many analyses

- Uncertainty is weight for analysis



# Slide 3 

Results - Lots of variation






# Figure Creators

```r
library(ggplot2)
library(dplyr)

xseq <- 42:52

seq_1a <- data.frame(
    entry = c("", "A", "T", "T", "A", "T", "G", "A", "", "", ""),
    score = c(1, 0.99, 0.999, 0.9999, 0.9998,  0.55, 0.97, 0.96, 1, 1, 1),
    x = xseq,
    label = "1a"
)


seq_1b <- data.frame(
    entry = c("", "A", "T", "T", "A", "C", "G", "A", "", "", ""),
    score = c(1, 0.95, 0.99, 0.99968, 1,  0.65, 0.99, 0.99, 1, 1, 1),
    x = xseq,
    label = "1b"
)


seq_2 <- data.frame(
    entry = c("", "", "", "T", "A",   "T",  "G",  "A", "T",  "C", ""),
    score = c(1,   1,  1, 0.9,  0.999,  0.995, 0.99, 0.99,   0.9992, 0.99,  1),
    x = xseq,
    label = "2"
)


seq_3 <- data.frame(
    entry = c("T", "A",  "T",  "T",  "A",   "C",  "G",  "A", "T",  "C", ""),
    score = c(1,   0.6, 0.55, 0.97, 0.92,  0.51, 0.99, 0.99,   1, 0.99,  1),
    x = xseq,
    label = "3"
)


seq_4 <- data.frame(
    entry = c("", "", "A", "T", "T",   "C",  "G",  "A", "T",  "C", "A"),
    score = c(1,   1,  1,  0.9,   1,  0.995, 0.99, 0.99,   0.99975, 0.99,  1),
    x = xseq,
    label = "4"
)

seq_con <- data.frame(
    entry = c("N", "A", "T", "T", "A", "C", "G", "A", "T", "C", "N"),
    score = rep(1, length(xseq)),
    x = xseq,
    label = "Conseq"
)


these_seqs <- bind_rows(
    seq_1a, seq_1b, seq_2, seq_3, seq_4, seq_con
)
these_seqs$label = factor(these_seqs$label, 
    ordered = TRUE, 
    levels = c("Conseq", "4", "3", "2", "1b", "1a"))

 ggplot(these_seqs) + 
        aes(x = x, y = label, label = entry, fill = score) +  
        theme_minimal() +  
        theme(plot.title = element_text(hjust = 0.5),  
            legend.position = "bottom", 
            text = element_text(size = 21), 
            legend.text = element_text(size = 11)) +  
        geom_tile() + 
        geom_text(size = 10) +  
        scale_fill_gradient(low = "red", high = "white") + 
        scale_x_continuous(breaks = 42:52) + 
        labs(y = "      Read Number", 
            x = "Nucleotide Location",  
            fill = "Quality", 
            title = "SAM file") + 
        geom_hline(yintercept = 1.5, size = 1)  
```

```r
# Quality scores
q <- data.frame(Q = (0:40), err_prob = 10^(-(0:40)/10)) %>%
    mutate(qual = 1 - err_prob)
q$phred <- sapply(q$Q + 33, intToUtf8)

these_seqs$phred <- sapply(these_seqs$score, 
    function(x) {
        q$phred[which.min(abs(q$qual - x))]
})

sapply(unique(these_seqs$label), 
    function(x) paste0(these_seqs$phred[these_seqs$label == x & these_seqs$entry != ""], collapse = ""))

sapply(unique(these_seqs$label), 
    function(x) paste0(these_seqs$entry[these_seqs$label == x & these_seqs$entry != ""], collapse = ""))

```

```r
# Quality Score Plot
library(stringr)
ggplot(q, aes(x = Q, y = qual)) +  
    theme_minimal() + 
    theme(text = element_text(size = 21)) + 
    geom_point(size = 3) + 
    scale_x_continuous(breaks = 0:40, 
        labels = q$phred, minor_breaks = NULL) +
    labs(x = NULL, y = "Quality")
  
```

```r
# Uncertainty Matrix
these_seqs2 <- filter(these_seqs, label != "Conseq")
seq_mat <- matrix(0, ncol = 11, nrow = 4)
rownames(seq_mat) <- c("A", "T", "C", "G")
colnames(seq_mat) <- xseq

for (i in seq_len(nrow(these_seqs2))) {
    score = these_seqs2$score[i]
    if (these_seqs2$label[i] %in% c("1a", "1b")) score <- score/2

    col <- which(colnames(seq_mat) == these_seqs2$x[i])
    row <- which(rownames(seq_mat) == these_seqs2$entry[i])
    seq_mat[row, col] <- seq_mat[row, col] + score
}

```



