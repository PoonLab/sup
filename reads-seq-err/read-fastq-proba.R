#
# Read FASTQ files, extract and translate error probabilities
# from Phred Q-scores
#
library(tidyr)
library(dplyr)
library(ggplot2) ; theme_set(theme_bw())
library(gridExtra)

set.seed(1234)


# ==== FUNCTIONS ====

# Converts the ASCII score to an integer Q:
#(http://drive5.com/usearch/manual/quality_score.html)
translate_to_Q <- function(x) {
    return(utf8ToInt(x)-33)
}

# Converts Q to the error probability:
# https://en.wikipedia.org/wiki/Phred_quality_score
translate_to_P <- function(x){
    Q <- translate_to_Q(x)
    return(10^(-Q/10))
}


fastq_seq_proba <- function(filename) {  # filename = 'novaseq_reads_R1.fastq'

    x <- readLines(filename)

    # Identify the lines of seq ID
    q.seqid <- which(grepl("^@",x))

    seqid <- x[q.seqid]
    seq   <- x[q.seqid+1]
    phred <- x[q.seqid+3]

    # Vector format (more convenient):
    phred <- strsplit(x = phred, split='')

    # Convert all reads to Q scores and probabilities:
    Q <- list()
    P <- list()
    for(i in 1:length(phred)){
        Q[[i]] <- sapply(phred[[i]], translate_to_Q)
        P[[i]] <- sapply(phred[[i]], translate_to_P)
    }

    return(list(seqid     = seqid,
                seq       = seq,
                err.proba = P,
                Q.score   = Q))

}



#' Read the quality FASTQ file and
#' translate ASCII scores to Q scores and
#' error probabilities
#' @param filename String. FASTQ file name.
#'
fastq_to_QP <- function(filename) {

    x <- readLines(filename)

    # Identify the line of quality scores
    # (comes after the '+' separator):
    q.lines <- which(x=='+')+1

    # Keep the scores only:
    q.reads <- x[q.lines]

    # Vector format (more convenient):
    qv <- strsplit(x = q.reads, split='')


    # Convert all reads to Q scores and probabilities:
    Q <- list()
    P <- list()
    for(i in 1:length(qv)){
        Q[[i]] <- sapply(qv[[i]], translate_to_Q)
        P[[i]] <- sapply(qv[[i]], translate_to_P)
    }

    allQ <- do.call('c',Q)
    allP <- do.call('c',P)

    return(list(Q = Q,
                P = P,
                allQ = allQ,
                allP = allP))
}


#' Plot the Q-scores for the ith read.
plot_Q_read <- function(i, Q) {
    qi = Q[[i]]
    plot(x=1:length(qi), y=qi,
         xlab = 'position', ylab = 'Q-score',
         ylim=c(0,45),
         main = paste('Read',i),
         las = 1,
         pch=16, cex = 2,
         col=rgb(0,0,0,0.5))
    grid()
}


plot_P_read <- function(i, P) {
    y = P[[i]]
    plot(x=1:length(y), y=y,
         xlab = 'position',
         ylab = 'Error probability',
         ylim=c(10^-5,1),
         main = paste('Read',i),
         las = 1,
         log = 'y',
         pch=16, cex = 1,
         col=rgb(0,0,0,0.5))
    grid()
}


reads_to_df <- function(res) {
    P <- res$P
    tmp <- list()
    for(i in 1:length(P)){
        tmp[[i]] <- data.frame(pos = 1:length(P[[i]]),
                               p = P[[i]],
                               read = i)
    }
    df <- do.call('rbind',tmp)
    return(df)
}


plot_P_pos <- function(res, title, ci = 0.99) {

    df <- reads_to_df(res)
    qt <- c(0.5-ci/2, 0.5+ci/2 )

    dfs <- df %>%
        group_by(pos) %>%
        summarise(p.mean = mean(p),
                  p.lo = quantile(p, probs = qt[1]),
                  p.hi = quantile(p, probs = qt[2])
        )

    g <- dfs %>%
        ggplot(aes(x=pos)) +
        # geom_point(aes(y=p.mean))+
        geom_ribbon(aes(ymin=p.lo, ymax=p.hi),
                    fill = 'blue',
                    alpha=0.3)+
        geom_line(aes(y=p.mean), alpha=0.7)+
        scale_y_log10(breaks = 10^c(-4:0), limits=10^c(-4,0))+
        xlab('position')+
        ylab('Error Probability')+
        ggtitle(title, paste('Mean, CI =',ci*100,'%'))

    return(g)
}

combo <- function(f) {
    z <- f %>%
        fastq_to_QP() %>%
        plot_P_pos(title = f)
    return(z)
}

# ==== RUN ====

z <- fastq_seq_proba(filename = 'novaseq_reads_R1.fastq')


i = 7
z$seqid[i]
z$seq[i]
z$err.proba[i]


# Read the FASTQ file:
fnames <- system('ls *.fastq', intern = TRUE)

# Translate to error probabilities
# and plot them by position:
plist <- lapply(fnames, combo)

pdf('plot-errorProba.pdf', width=12, height = 12)
do.call("grid.arrange", c(plist, ncol=2))
dev.off()

