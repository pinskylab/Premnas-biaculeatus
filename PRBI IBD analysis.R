#Premnas biaculeatus Genotypes and Analysis

library(genepop)

#Read in data

prbi_genalex = read.csv("PRBI_GenAlEx_2009-11-12.csv")
prbi_genepop = read.csv("PRBI_2009-09-07.txt")

#Calculate pairwise Fst between populations
#pairs=TRUE to get a pairwise Fst matrix

Fst("PRBI_2009-09-07.txt", pairs = TRUE, outputFile = "prbi_genepop_fst.txt", 
    dataType = "Diploid", verbose = interactive())

#Resulting Fst matrices found in "prbi_genepop_fst.txt" and "prbi_genepop_fst.txt.MIG" 

#Calculate geographic distances between populations using geodist package

install.packages("geodist")
library(geodist)

#Create vector for longitude and latitude of populations
#Listed in sequential order of populations from pop 1-pop 22

prbi_long <- c(124.2489589, 124.177941, 123.3179946, 123.4794616, 
               123.5785325, 123.6545571, 123.8037213, 124.0285911,
               124.0109535, 124.032865, 125.0251015, 124.6560577)
prbi_long

prbi_lat <- c(10.22834459, 10.17666235, 9.414770978, 9.616869822,
              9.848377696, 10.07330622, 10.22512305, 10.48777602,
              10.76495538, 10.9440459, 10.01649883, 10.54725344)
prbi_lat


#Use geodist_vec to calculate pairwise distances
#Outputs matrix of geodesic distances in meters- labeled as prbi_geodistance

prbi_geodistance_meters <- geodist_vec(prbi_long, prbi_lat, paired=FALSE, sequential=FALSE,
            pad=FALSE, measure="geodesic")
prbi_geodistance_meters 

#Convert geographic distance from meters to kilometers

prbi_geodistance_km <- prbi_geodistance_meters/1000
prbi_geodistance_km

#Create a matrix for pairwise Fst

prbi_fst = read.csv("Pairwise Fst Matrix.csv")
prbi_fst_matrix <- as.matrix(prbi_fst)

#Plot pairwise Fst matrix vs. geographic distance matrix

plot(prbi_geodistance_km, prbi_fst_matrix, ylim=c(-0.25, 0.25), 
     xlab="Pairwise Geographic Distance (km)", ylab="Pairwise Fst", 
     main="Pairwise Fst vs. Geographic Distance", pch=20)

#IBD analysis using Genepop function
library(genepop)

#Tried this but R aborted and terminated
#Not sure on the format of the input file
ibd("prbi_genepop_ibdinput.txt", outputFile="prbi_genepop_ibdanalysis.txt", 
    dataType="Diploid", statistic="F/(1-F)", geographicScale="1D", CIcoverage=0.95,
    testPoint=0, mantelPermutations=1000, verbose = interactive())


