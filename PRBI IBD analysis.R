#Premnas biaculeatus Genotypes and Analysis

library(genepop)

##Read in data

prbi_genalex = read.csv("PRBI_GenAlEx_2009-11-12.csv")
prbi_genepop = read.csv("PRBI_2009-09-07.txt")

##Test for Hardy-Weinberg Equilibrium

test_HW("PRBI_2009-09-07.txt", outputFile="prbi_hw.txt")

#HW test results in file "prbi_hw.txt"


##Calculate pairwise Fst between populations
#pairs=TRUE to get a pairwise Fst matrix

Fst("PRBI_2009-09-07.txt", pairs = TRUE, outputFile = "prbi_genepop_fst.csv", 
    dataType = "Diploid", verbose = interactive())

#Resulting Fst matrices found in "prbi_genepop_fst.csv" and "prbi_genepop_fst.csv.MIG" 


##Calculate geographic distances between populations using geodist package

install.packages("geodist")
library(geodist)
install.packages("tidyverse")
library(tidyverse)

#Import avg. lat and long from csv file

prbi_avglatlong <- read.csv("PRBI_AvgLatLong.csv")

#Rename and reorder columns
prbi_avglatlong <- rename(prbi_avglatlong, latitude = Avg..Latitude, 
                                longitude = Avg..Longitude)

prbi_avglatlong <- prbi_avglatlong[c("longitude", "latitude", "Population")]

prbi_avglatlong_matrix <- as.matrix(prbi_avglatlong)

#Use geodist function to calculate distances (output is in meters)
#Getting error- says it can't determine longitude and latitude columns
prbi_distance_meters <- geodist(prbi_avglatlong, measure="geodesic")

#Use vector version of Geodist instead
#Create vector for longitude and latitude of populations from prbi_avglatlong
#Listed in sequential order of populations from pop 1-pop 22

prbi_long <- prbi_avglatlong$longitude
prbi_long

prbi_lat <- prbi_avglatlong$latitude
prbi_lat

#Use geodist_vec to calculate pairwise distances
#Outputs matrix of geodesic distances in meters- labeled as prbi_geodistance_meters

prbi_geodistance_meters <- geodist_vec(prbi_long, prbi_lat, paired=FALSE, sequential=FALSE,
            pad=FALSE, measure="geodesic")
prbi_geodistance_meters 

#Convert geographic distance from meters to kilometers

prbi_geodistance_km <- prbi_geodistance_meters/1000
prbi_geodistance_km

#Relabel matrix columns and rows to be Populations 1-22
rownames(prbi_geodistance_km) <- c(1, 2, 7, 8, 9, 10, 11, 13, 14, 15, 19, 22)
colnames(prbi_geodistance_km) <- c(1, 2, 7, 8, 9, 10, 11, 13, 14, 15, 19, 22)
prbi_geodistance_km

#Remove diagonals/upper triangle of the geographic distance matrix
prbi_geodistance_km[upper.tri(prbi_geodistance_km, diag=T)] = NA



##Fst of all populations and all loci
#Trying different ways of loading in the fst matrix to have less formatting issues
#None worked so far
prbi_txt_fst <- read.table("prbi_genepop_fst.txt.MIG", header=T, skip=3, sep="")

library(data.table)
prbi_txt_fst <- fread("prbi_genepop_fst.txt.MIG")

prbi_txt_fst <- read.table("prbi_genepop_fst.csv", skip=270, sep="")

#Load in Fst matrix and edit formatting to match prbi_geodistance_km
prbi_fst <- read.csv("prbi_genepop_fst.csv.MIG")
prbi_fst <- as.matrix(prbi_fst)

#Delete unneeded rows 1, 2, and 14
prbi_fst <- prbi_fst[c(-1, -2, -14)]
#still have the issue of only one column holding all the Fst values

#Instead will import a csv file that has the Fst values in appropriate columns
prbi_fst <- read.csv("Pairwise Fst Matrix.csv")
prbi_fst <- as.matrix(prbi_fst)
prbi_fst <- prbi_fst[,-1]

colnames(prbi_fst) <- c(1, 2, 7, 8, 9, 10, 11, 13, 14, 15, 19, 22)
rownames(prbi_fst) <- c(1, 2, 7, 8, 9, 10, 11, 13, 14, 15, 19, 22)

##Linearize Fst and Plot Linear Fst vs. Geographic Distance

#Linearize Fst (Fst/(1-Fst))
prbi_fstlin = as.matrix(prbi_fst/(1-prbi_fst))

#Plot Linear Fst vs. Geographic Distance
plot(prbi_geodistance_km, prbi_fstlin, ylim=c(-0.25, 0.25), 
     xlab="Pairwise Geographic Distance (km)", ylab="Fst/(1-Fst)", 
     main="Linearized Fst vs. Geographic Distance", pch=20)


##Excluding ACHB9

##Calculate pairwise Fst between populations- using Genepop file excluding ACHB9
#pairs=TRUE to get a pairwise Fst matrix

Fst("PRBI_2009-09-07_noACH_B9.txt", pairs = TRUE, outputFile = "prbi_genepop_noACHB9_fst.csv", 
    dataType = "Diploid", verbose = interactive())

#Resulting Fst matrices found in "prbi_genepop_noACHB9_fst.csv" and "prbi_genepop_noACHB9_fst.csv.MIG" 

#Create matrix for Fst excluding ACHB9
prbi_noACHB9_fst <- read.csv("All Pop No ACHB9 Fst Matrix.csv")
prbi_noACHB9_fst <- as.matrix(prbi_noACHB9_fst)
prbi_noACHB9_fst <- prbi_noACHB9_fst[,-1]

colnames(prbi_noACHB9_fst) <- c(1, 2, 7, 8, 9, 10, 11, 13, 14, 15, 19, 22)
rownames(prbi_noACHB9_fst) <- c(1, 2, 7, 8, 9, 10, 11, 13, 14, 15, 19, 22)

##Linearize Fst and Plot Linear Fst vs. Geographic Distance

#Linearize Fst (Fst/(1-Fst))
prbi_noACHB9_fstlin = as.matrix(prbi_noACHB9_fst/(1-prbi_noACHB9_fst))

#Plot Linear Fst vs. Geographic Distance
plot(prbi_geodistance_km, prbi_noACHB9_fstlin, 
     xlab="Pairwise Geographic Distance (km)", ylab="Fst/(1-Fst)", 
     main="All Pop No ACHB9 Linearized Fst vs. Geographic Distance", pch=20)



##Extra/unused:
##Create a matrix for pairwise Fst

prbi_fst = read.csv("Pairwise Fst Matrix.csv")
prbi_fst_matrix <- as.matrix(prbi_fst)

##Plot pairwise Fst matrix vs. geographic distance matrix

plot(prbi_geodistance_km, prbi_fst_matrix, ylim=c(-0.25, 0.25), 
     xlab="Pairwise Geographic Distance (km)", ylab="Pairwise Fst", 
     main="Pairwise Fst vs. Geographic Distance", pch=20)


##IBD analysis using Genepop function
library(genepop)

#Tried this but R aborted and terminated
#Not sure on the format of the input file
ibd("prbi_genepop_ibdinput.txt", outputFile="prbi_genepop_ibdanalysis.txt", 
    dataType="Diploid", statistic="F/(1-F)", geographicScale="1D", CIcoverage=0.95,
    testPoint=0, mantelPermutations=1000, verbose = interactive())

#Trying a version of the ibd function that allows a geo distance matrix
#Error- unused argument geoDistFile- didn't take it as an input
ibd("PRBI_2009-09-07.txt", geoDistFile="prbi_genepop_geomatrix.txt", 
    outputFile="prbi_genepop_ibdanalysis.txt", dataType="Diploid", 
    statistic="F/(1-F)", geographicScale="1D", CIcoverage=0.95,
    testPoint=0, mantelPermutations=1000, verbose = interactive())





