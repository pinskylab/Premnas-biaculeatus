#Premnas biaculeatus Genotypes and Analysis

library(genepop)

#Read in data

prbi_genalex = read.csv("PRBI_GenAlEx_2009-11-12.csv")
prbi_genepop = read.csv("PRBI_2009-09-07.txt")

#Calculate pairwise Fst between populations
#Unsure on pairs= and sizes= arguments- what does True and False mean for each?

Fst("PRBI_2009-09-07.txt", pairs = TRUE, outputFile = "prbi_genepop_fst.txt", 
    dataType = "Diploid", verbose = interactive())

#Resulting Fst matrices found in "prbi_genepop_fst.txt" and "prbi_genepop_fst.txt.MIG" 