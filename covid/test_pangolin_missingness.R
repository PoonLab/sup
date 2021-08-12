library(stringr)

sampled <- list.files("data/sampled_covid", pattern = "*.fasta")

# Import one sequence

for(thisfile in seq_along(sampled)){
    testname <- sampled[thisfile]
    testfasta <- readLines(con = paste0("data/sampled_covid/", testname),
        n = 2)[2]

    alph <- c("A", "T", "C", "G")

    atcg <- str_count(testfasta, alph)
    n <- str_count(testfasta, "N")

    atcg_rates <- atcg/(nchar(testfasta) - n)
    barplot(atcg_rates, names.arg = alph)

    # Step 1: Remove Ns with random letters (unrealistic) ----
    test2 <- testfasta
    counter <- 0
    while(str_count(test2, "N") > 0){
        counter <- counter + 1
        if(counter > nchar(test2)) break

        newletter <- sample(alph, 1, prob = atcg_rates)
        test2 <- sub(pattern = "N", replacement = newletter, test2)
    }
    str_count(test2, "N") # success




    # Step 2: Randomly add N's ----

    # Set up number of missing values
    Nvals <- round(nchar(test2) * c(0.001, 0.005, 0.01, 0.05, 0.1, 0.25, 0.5))
    # Number of samples with each amount of missingness
    Nsamps <- 1000

    # Prep the data frame
    samps <- expand.grid(id = 1:Nsamps, NVals = Nvals)
    samps$seq <- test2

    # Replace characters in data frame with Ns
    for(i in 1:nrow(samps)){
        cat("\r", i)
        newns <- sample(1:nchar(samps$seq[i]), samps$NVals[i], replace = FALSE)
        thisseq <- strsplit(samps$seq[i], "")[[1]]
        thisseq[newns] <- "N"
        samps$seq[i] <- paste0(thisseq, collapse = "")

    }; cat("\n")

    all.equal(str_count(samps$seq, "N"), samps$NVals)

    refseq <- paste0("> 0.1\n", test2, "\n", collapse = "")

    allseq <- paste0(refseq,
        paste0("> ", samps$NVals, ".", samps$id, "\n", samps$seq, "\n", collapse = ""),
        collapse = ""
    )

    writeLines(allseq,
        paste0("/mnt/BCC20BCCC20B8A3A/Fasta Files/sra-downloads/files/", testname))

    testcsv <- paste0(strsplit(testname, "\\_")[[1]][1], ".csv", collapse = "")
    print(paste0("pangolin ", testname, " --outfile pangolin_results/", testcsv, collapse = ""))
}

# conda activate pangolin
# pangolin missingness2.fasta --outdir pangolin_results


if(FALSE){ # testing the two extremes
    library(stringr)
    sampled <- c("ERR4693034_sampled.fasta", "ERR4693061_sampled.fasta")
    # Import one sequence

    for(thisfile in seq_along(sampled)){
        testname <- sampled[thisfile]
        testfasta <- readLines(con = paste0("data/sampled_covid/", testname),
            n = 2)[2]

        alph <- c("A", "T", "C", "G")

        atcg <- str_count(testfasta, alph)
        n <- str_count(testfasta, "N")

        atcg_rates <- atcg/(nchar(testfasta) - n)
        #barplot(atcg_rates, names.arg = alph)

        # Step 1: Remove Ns with random letters (unrealistic) ----
        test2 <- testfasta
        counter <- 0
        while(str_count(test2, "N") > 0){
            counter <- counter + 1
            if(counter > nchar(test2)) break

            newletter <- sample(alph, 1, prob = atcg_rates)
            test2 <- sub(pattern = "N", replacement = newletter, test2)
        }
        str_count(test2, "N") # success




        # Step 2: Randomly add N's ----

        # Set up number of missing values
        Nvals <- round(nchar(test2) * c(0.005, 0.01, 0.05,  
            seq(0.05, 0.2, length.out = 40), 
            0.25, 0.3, 0.4, 0.45, 0.5))
        # Number of samples with each amount of missingness
        Nsamps <- 100

        # Prep the data frame
        samps <- expand.grid(id = 1:Nsamps, NVals = Nvals)
        samps$seq <- test2

        # Replace characters in data frame with Ns
        for(i in 1:nrow(samps)){
            cat("\r", i, "/", nrow(samps))
            newns <- sample(1:nchar(samps$seq[i]), samps$NVals[i],
                replace = FALSE)
            thisseq <- strsplit(samps$seq[i], "")[[1]]
            thisseq[newns] <- "N"
            samps$seq[i] <- paste0(thisseq, collapse = "")

        }; cat("\n")

        #all.equal(str_count(samps$seq, "N"), samps$NVals)

        refseq <- paste0("> 0.1\n", test2, "\n", collapse = "")

        allseq <- paste0(refseq,
            paste0("> ", samps$NVals, ".", samps$id, "\n", samps$seq, "\n",
                collapse = ""),
            collapse = ""
        )

        thisstripped <- paste0(strsplit(testname, "\\_")[[1]][1],
            "-2", collapse = "")

        newfasta <- paste0("/mnt/BCC20BCCC20B8A3A/Fasta Files/sra-downloads/files/",
            thisstripped, ".fasta")
        print(newfasta)


        testcsv <- paste0(thisstripped, ".csv", collapse = "")
        print(paste0("pangolin ", newfasta, " --outfile pangolin_results/",
            testcsv, collapse = ""))

        writeLines(allseq, con = newfasta)
    }
}


