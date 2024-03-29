library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(here)
library(patchwork)
library(scales)
library(xtable)
source(here("covid", "aux_funk.R"))



#```{r load_sampled, include=FALSE}
# Read in CSV files
lins <- read.csv(here("data", "output", "lins.csv"),
    header = TRUE, stringsAsFactors = FALSE)

taxons <- sort(unique(lins$taxon))
#```

#```{r load_ordered, include=FALSE, eval=FALSE}
ords <- lapply(
    list.files(here("data/pangordlineages"),
        full.names = TRUE),
    read.csv
    ) %>%
    bind_rows() %>%
    separate(col = taxon,
        into = c("dummy", "acc", "lik", "subs"),
        sep = "_") %>%
    mutate(nsubs = str_count(subs, "2"),
        lik = as.numeric(lik)) %>%
    group_by(acc) %>%
    mutate(order = 1:n()) %>%
    ungroup()

#```




#```{r sampled_bars, fig.height=16, fig.width=13, fig.cap='\\label{fig:covidcalls}Distribution of called lineages from pangolin. Red bars indicate the lineage of the most probable sequence. Any lineage with fewer than 100 observations was grouped into the "Other" category.'}
max_label <- 250
other_label <- 100

pdf(here("ms", "figs", "sampled_bars.pdf"), width = 14, height = 17)
par(mfrow = c(45, 1), mar = c(0.05, 7.75, 0.05, 0.05))
if (exists("seq_info")) rm(seq_info)
for (i in seq_along(taxons)) {
    pang <- lins[lins$taxon == taxons[i], ]
    called <- pang$lineage[pang$sample == 0][1]
    pangtab <- sort(table(pang$lineage), decreasing = TRUE)

    # Prep the data for a nicely formatted table
    # Subtract one because of the conseq.
    seq_info_i <- data.frame(
        called = called,
        mode = names(pangtab)[1],
            mode_n = pangtab[1] - 1,
            perc = round(100 * (pangtab[1] - 1) / (sum(pangtab) - 1), 2),
        runner_up = names(pangtab)[2],
            ru_n = pangtab[2],
        unique = length(pangtab), atoms = sum(pangtab == 1))
    seq_info_i$taxon <- taxons[i]

    if (!exists("seq_info")) {
        seq_info <- seq_info_i
    } else {
        seq_info <- bind_rows(seq_info, seq_info_i)
    }

    colvec <- rep("grey", length(pangtab))
    colvec[which(names(pangtab) == called)] <- "red"

    n <- sum(pangtab > max_label)
    if (n > 1) {
        add_other <- FALSE
        if (sum(pangtab < other_label) > 10) {
            add_other <- TRUE
            other_count <- sum(pangtab <= other_label)
            pangtab <- c(pangtab[pangtab > other_label],
                c("other" = sum(pangtab[pangtab <= other_label])))
            colvec[which(names(pangtab) == "other")] <- "black"
        }
        barlabx <- c(0, cumsum(pangtab[1:(n - 1)])) +
            pangtab[1:n] / 2
        barlabels <- names(pangtab)[1:n]
        barlens <- sapply(gregexpr("\\.", barlabels), length)
        for (j in seq_along(barlabels)) {
            if (pangtab[j] < 400 & barlens[j] >= 2) {
                barsplit <- strsplit(barlabels[j], split = "\\.")[[1]]
                barn <- length(barsplit)
                half <- floor(barn / 2)
                barlabels[j] <- paste0(
                    paste(barsplit[1:half], collapse = "."),
                    ".\n",
                    paste(barsplit[(half + 1):barn], collapse = ".")
                )
            }
        }

        barplot(as.matrix(pangtab),
            col = colvec, hori = TRUE, axes = FALSE)
        text(barlabx, 0.7, barlabels, cex = 1.5)
        if (add_other) {
            text(x = sum(pangtab) - pangtab["other"] / 2,
            y = 0.7, col = "white", cex = 1.5,
            label = paste0("Others:\n", other_count))
        }
        mtext(side = 2, cex = 1, las = 1,
            text = paste(substr(taxons[i], 1, 3),
                substr(taxons[i], 4, 20), sep = ""))
        abline(v = seq(0, 10000, 1000), lty = 2)
        "pretty_labels <- seq(0, sum(pangtab),
                by = ifelse(sum(pangtab) < 2000, 100, 1000))
        mtext(side = 1,
            at = pretty_labels,
            text = pretty_labels,
            line = 0,
            cex = 0.75
        )"
    }
}
dev.off()
#```


# \dc{I'd like to suggest an alternative to Figure 1. A histogram where the x-axis is the proportion of resamples that are assigned the same lineage as the consensus sequence. This way we would have a figure that represent \emph{all} the data (not truncated to those that have more than 250 observations); this may also be more to the point: taking sequencing uncertainty into account, some (many?) sequences cannot be so confidently assigned to a single consensus sequence. We could even add another panel (so that would be a 2-panel figure) with a histogram again, but this time the x-axis shows the number of lineages associated to the sequences with probability, say, greater than 5\%.}
prop_correct <- prop_second <- prop_third <- double(length(taxons))
atoms <- len_unique <- one_v_two <- double(length(taxons))

for (i in seq_along(taxons)) {
    pang <- lins[lins$taxon == taxons[i], ]
    called <- pang$lineage[pang$sample == 0][1]
    pangtab <- sort(table(pang$lineage), decreasing = TRUE)

    prop_correct[i] <- pangtab[names(pangtab) == called]/sum(pangtab)
    prop_second[i] <- pangtab[2] / sum(pangtab)
    prop_third[i] <- pangtab[3] / sum(pangtab)

    atoms[i] <- sum(pangtab == 1)
    len_unique[i] <- length(pangtab)
    one_v_two[i] <- pangtab[1] / pangtab[2]
}

bins = seq(0, 1, 0.01)
pdf(here("ms", "figs", "histograms.pdf"), width = 9, height = 6)
par(mfrow = c(2, 3))
hist(prop_correct, breaks = bins,
    main = "Proportion of resamples agreeing\nwith consensus sequence,\nnot including 0s",
    xlab = "Proportion")
hist(prop_second[prop_second > 0], breaks = bins, 
    main = "Proportion of resamples agreeing\nwith second most common lineage\nnot including 0s",
    xlab = "Proportion")
hist(prop_third[prop_third > 0], breaks = bins,
    main = "Proportion of resamples agreeing\nwith third most common lineage\nnot including 0s",
    xlab = "Proportion")
hist(atoms, breaks = 40,
    main = "Number of atoms sampled",
    xlab = "Count")
hist(len_unique, breaks = 40,
    main = "Number of unique lineages sampled",
    xlab = "Count")
hist(one_v_two, breaks = 40,
    main = "Ratio of most common versus\nsecond most common lineage")
par(mfrow = c(1,1))
dev.off()

hinset <- ggplot(mapping = aes(x = len_unique)) + 
    theme_bw() +
    theme(axis.text = element_text(size = 6), axis.title = element_text(size = 8)) + 
    geom_histogram(colour = "black", fill = "lightgrey") +
    labs(x = "Number of unique lineages in \neach set of resampled sequences",
        y = "Count")

pdf(here("ms", "figs", "prop_correct.pdf"), width = 6, height = 4)
ggplot(mapping = aes(x = prop_correct[prop_correct > 0.2 & prop_correct < 1])) + 
    geom_histogram(binwidth=0.025, center = 0.0125,
        colour = "black", fill = "lightgrey") +
    theme_bw() +
    labs(x = "Proportion of resampled sequences assigned to the same lineage \nas the consensus sequence, not including <1% or exactly 100%",
        y = "Count") +
    annotation_custom(ggplotGrob(hinset),
        xmin = 0.3, xmax = 0.7,
        ymin = 15, ymax = 65)
dev.off()

hist(prop_correct[prop_correct > 0 & prop_correct < 1])

#```{r RTT_slope, warning=FALSE, results='hide', fig.cap="\\label{fig:RTT_slope}Clock rates (slope) and 95% Confidence Intervals for the collections of re-sampled sequences. The red line and red shaded region are the clock rate and 95% CI for the consensus sequences. The purple points and error bars are the clock rates and error intervals (either Bayesian Credible Interval or Highest Posterior Probability) from published studies, as labelled. The re-sampled sequences are in line with the consensus sequences as well as the published sequences, but represent a much larger variation due to the uncertainty in the original genome sequences."}


# raw (conseqs), tree_dirs (resampled), and 
    # lit_clock (published)
load(here("RTT", "RTT_results.RData"))
lit_clock_x <- (max(as.numeric(tree_dirs$rank))+1):
                (max(as.numeric(tree_dirs$rank)) + 
                    nrow(lit_clock))

ggplot() +
    theme_bw() +
    geom_hline(yintercept = 0, colour = "darkgrey") +
    geom_hline(yintercept = raw$slope, col = "red") +
    geom_rect(
        mapping = aes(xmin = -10, xmax = 100,
            ymin = raw$slope - 1.96*raw$sd, ymax = raw$slope + 1.96*raw$sd),
        fill = "red", alpha = 0.2) +
    geom_errorbar(
        mapping = aes(x = as.numeric(rank),
            ymin = slope - 1.96*sd,
            ymax = slope + 1.96*sd),
        data = tree_dirs) +
    geom_point(
        mapping = aes(x = as.numeric(rank), y = slope),
        data = tree_dirs) +
    geom_errorbar(
        mapping = aes(x = lit_clock_x,
            ymin = as.numeric(lo95),
            ymax = as.numeric(hi95)),
        data = lit_clock,
        colour = "darkorchid") +
    geom_point(
        data = lit_clock,
        mapping = aes(x = lit_clock_x, 
            y = as.numeric(clock)),
        colour = "darkorchid"
        ) +
    geom_text(data = lit_clock,
        mapping = aes(
            y = as.numeric(hi95), 
            x = lit_clock_x,
            label = study),
        colour = "black", 
        angle = 90,
        hjust = -0.05,
        size = 2.5
        ) +
    scale_colour_viridis_d() +
    coord_cartesian(xlim = c(0, 50), ylim = c(-0.001, 0.0035)) +
    labs(x = "Index (Ordered by Slope)", 
        y = "Slope +/- 1.96 SD",
        title = "Clock rate estimates from collections of re-sampled genomes",
        subtitle = "Compared to conseqs (red) and published estimates (purple)") +
    theme(legend.position = "none") +
    scale_y_continuous(
        breaks = pretty(seq(-0.001, 0.0035, by = 0.00001), 
            n = 15)
        )

if (FALSE) {
    raw
    mean(tree_dirs$slope)
    sd(tree_dirs$slope)
    qqnorm(tree_dirs$slope)
    qqline(tree_dirs$slope)
}
#```



#```{r big_kable, fig.cap='\\label{tab:pango}Results for re-sampling analysis all accession numbers in our analysis. The "Conseq" column refers to the lineage designated to the consensus sequence. "Mode" refers to the most common lineage across all re-samples, with "Mode #" indicating the number of re-samples with this lineage. "Runner Up" and "RU #" are the second most likely lineage designation and the number with that designation, respectively. "Unique" is the total number of unique lineage designations, and "Atoms" are the total number of lineage designations that were only observed for a single re-sample.'}
if (file.exists(here("data", "pangoLEARN_lineage-recall-report.csv"))) {
    lin_report <- read.csv(here("data", 
        "pangoLEARN_lineage-recall-report.csv"))
} else { # download from pagoLEARN
    lin_report <- read.csv("https://raw.githubusercontent.com/cov-lineages/pangoLEARN/master/pangoLEARN/data/lineage_recall_report.csv")
    write.csv(lin_report, file = here("data", "pangoLEARN_lineage-recall-report.csv"))
}

lin_report <- select(lin_report, lineage, F1 = f1_score) %>%
    mutate(F1 = round(F1*100, 2))


seq_info$taxon <- taxons
seq_info <- arrange(seq_info, mode, mode_n) %>%
    select(taxon, everything())
seq_info <- rename(seq_info, Accession = taxon, Conseq = called,
    Mode = mode, `Mode #` = mode_n, `Mode #/N` = perc,
    `Runner-Up` = runner_up, `RU #` = ru_n,
    Unique = unique, Atoms = atoms)

seq_info <- left_join(seq_info, lin_report, by = c("Conseq" = "lineage"))

knitr::kable(seq_info, row.names = FALSE)
#```



# Generate table of accession numbers for Resampling
resampled_acc <- unique(lins$taxon)
# Pad with NAs to avoid recycling
resampled_acc <- resampled_acc[1:(ceiling(length(resampled_acc)/6)*6)]
print(xtable(matrix(resampled_acc, ncol = 6)), include.rownames=FALSE)

# Generate table of accession numbers for RTT
rtt_acc <- here("RTT", "samDONE") %>%
    list.files() %>%
    parse_accession() %>%
    unique()
# Pad with NAs to avoid recycling
rtt_acc <- rtt_acc[1:(ceiling(length(rtt_acc)/6)*6)]
print(xtable(matrix(rtt_acc, ncol = 6)), include.rownames=FALSE)



# Data for application section
lindata <- lins %>%
    group_by(taxon) %>%
    summarise(agree = mean(lineage[sample == 0] == lineage[sample != 0]),
        max_lin = names(sort(table(lineage), decreasing = TRUE))[1],
        plurality = lineage[sample == 0],
        none = sum(lineage[sample != 0] == "None"),
        nuns = lineage[sample == 0] == "None") # Plurality none -> "None"s -> nuns

nrow(lindata)
mean(lindata$plurality == "B.1.1.7")
sum(lindata$max_lin == "B.1.1.7")
sum(lindata$agree == 0)
sum(lindata$agree == 1)
lindata[lindata$agree == 0, ]  %>% arrange(plurality) %>% print(n = Inf) 
lindata[lindata$agree == 1, ]  %>% arrange(plurality) %>% print(n = Inf) 
lindata[lindata$plurality == "B.1.1.7",]
arrange(lindata, agree) %>% print(n = 60)
arrange(lindata, -none) %>% print(n = 100)
lindata[lindata$nuns, ]


