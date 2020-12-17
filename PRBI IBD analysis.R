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

#Load in Fst matrix
prbi_fst <- read.csv("prbi_genepop_fst.csv.MIG")
prbi_fst <- as.matrix(prbi_fst)

##Linearize Fst and Plot Linear Fst vs. Geographic Distance

#Read in Fst and Geographic Distance Matrices
#These ones only have information below the diagonal
#Need to use as.matrix to move this to matrix form (otherwise geo distance is not in a matrix)

prbi_fst_trimatrix = read.csv("Pairwise Fst Triangle Matrix.csv")
prbi_geodistance_trimatrix = read.csv("Geographic Distance Triangle Matrix.csv")


#Linearize Fst (Fst/(1-Fst))
prbi_fstlin = as.matrix(prbi_fst_trimatrix/(1-prbi_fst_trimatrix))

#Plot Linear Fst vs. Geographic Distance
plot(prbi_geodistance_trimatrix, prbi_fstlin)






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





