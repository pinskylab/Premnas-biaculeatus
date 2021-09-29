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

#Add regression to graph
library(smatr)

line_all <- line.cis(y=c(prbi_fstlin), x=c(prbi_geodistance_km))
abline(a=line_all$coef[1], b=line_all$coef[2], col="red", lwd=1)

#Mantel test

library(vegan)
mantel(prbi_geodistance_km, prbi_fstlin) 
#r: 0.1563 p: 0.146

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

#Add regression to graph

line_noACHB9 <- line.cis(y=c(prbi_noACHB9_fstlin), x=c(prbi_geodistance_km))
abline(a=line_noACHB9$coef[1], b=line_noACHB9$coef[2], col="red", lwd=1)

#Mantel test

mantel(prbi_geodistance_km, prbi_noACHB9_fstlin) 
#r: 0.1546 p: 0.167


##Excluding Population 19

#Create geographic distance matrix that excludes population 19

#Exclude pop 19 from the avg. lat long dataframe
prbi_no19_avglatlong <- prbi_avglatlong[-c(11),]

#Create vector for longitude and latitude of populations from prbi_no19_avglatlong
#Listed in sequential order of populations from pop 1-pop 22
prbi_no19_long <- prbi_no19_avglatlong$longitude
prbi_no19_long

prbi_no19_lat <- prbi_no19_avglatlong$latitude
prbi_no19_lat

#Use geodist_vec to calculate pairwise distances
#Outputs matrix of geodesic distances in meters- labeled as prbi_no19_geodist_meters
library(geodist)

prbi_no19_geodist_meters <- geodist_vec(prbi_no19_long, prbi_no19_lat, paired=FALSE, sequential=FALSE,
                                       pad=FALSE, measure="geodesic")
prbi_no19_geodist_meters 

#Convert geographic distance from meters to kilometers

prbi_no19_geodist_km <- prbi_no19_geodist_meters/1000
prbi_no19_geodist_km

#Relabel matrix columns and rows to be Populations 1-22
rownames(prbi_no19_geodist_km) <- c(1, 2, 7, 8, 9, 10, 11, 13, 14, 15, 22)
colnames(prbi_no19_geodist_km) <- c(1, 2, 7, 8, 9, 10, 11, 13, 14, 15, 22)
prbi_no19_geodist_km

#Remove diagonals/upper triangle of the geographic distance matrix
prbi_no19_geodist_km[upper.tri(prbi_no19_geodist_km, diag=T)] = NA


#Generate Fst matrices excluding population 19 and excluding ACHB9:

##Calculate pairwise Fst between populations- using Genepop file excluding pop 19 and ACHB9
#pairs=TRUE to get a pairwise Fst matrix
library(genepop)

Fst("PRBI_2009-09-07_no19_noACH_B9.txt", pairs = TRUE, outputFile = "prbi_genepop_no19_noACHB9_fst.csv", 
    dataType = "Diploid", verbose = interactive())

#Resulting Fst matrices found in "prbi_genepop_no19_noACHB9_fst.csv" and "prbi_genepop_no19_noACHB9_fst.csv.MIG" 

#Create matrix for Fst excluding ACHB9
prbi_no19_noACHB9_fst <- read.csv("No 19 No ACHB9 Fst Matrix.csv")
prbi_no19_noACHB9_fst <- as.matrix(prbi_no19_noACHB9_fst)
prbi_no19_noACHB9_fst <- prbi_no19_noACHB9_fst[,-1]

colnames(prbi_no19_noACHB9_fst) <- c(1, 2, 7, 8, 9, 10, 11, 13, 14, 15, 22)
rownames(prbi_no19_noACHB9_fst) <- c(1, 2, 7, 8, 9, 10, 11, 13, 14, 15, 22)

##Linearize Fst and Plot Linear Fst vs. Geographic Distance

#Linearize Fst (Fst/(1-Fst))
prbi_no19_noACHB9_fstlin = as.matrix(prbi_no19_noACHB9_fst/(1-prbi_no19_noACHB9_fst))

#Faster way to exclude pop 19 from the Fst matrix (delete the row and column where pop 19 data is):
prbi_no19_noACHB9_fstlin <- prbi_noACHB9_fstlin[-11, -11]

#Plot Linear Fst vs. Geographic Distance
plot(prbi_no19_geodist_km, prbi_no19_noACHB9_fstlin, 
     xlab="Pairwise Geographic Distance (km)", ylab="Fst/(1-Fst)", 
     main="No Pop 19 No ACHB9 Linearized Fst vs. Geographic Distance", pch=20)

#Add regression to graph

line_no19_noACHB9 <- line.cis(y=c(prbi_no19_noACHB9_fstlin), x=c(prbi_no19_geodist_km))
abline(a=line_no19_noACHB9$coef[1], b=line_no19_noACHB9$coef[2], col="red", lwd=1)

#Mantel test

mantel(prbi_no19_geodist_km, prbi_no19_noACHB9_fstlin) 
#r: 0.09557 p: 0.306


##Excluding population 19 but with all loci:
prbi_no19_fstlin <- prbi_fstlin[-11, -11]

#Plot Linear Fst vs. Geographic Distance
plot(prbi_no19_geodist_km, prbi_no19_fstlin, 
     xlab="Pairwise Geographic Distance (km)", ylab="Fst/(1-Fst)", 
     main="No Pop 19 Linearized Fst vs. Geographic Distance", pch=20)

#Add regression to graph

line_no19 <- line.cis(y=c(prbi_no19_fstlin), x=c(prbi_no19_geodist_km))
abline(a=line_no19$coef[1], b=line_no19$coef[2], col="red", lwd=1)

#Mantel test

mantel(prbi_no19_geodist_km, prbi_no19_fstlin) 
#r: 0.1158 p: 0.242



##Trial with a different Genepop input file found on GitHub:
#File "PRBI_genepop_2009-11-12.gen"

##Calculate pairwise Fst between populations
#pairs=TRUE to get a pairwise Fst matrix

Fst("PRBI_genepop_2009-11-12.gen.txt", pairs = TRUE, outputFile = "prbi_genepop_2009-11-12_fst.csv", 
    dataType = "Diploid", verbose = interactive())

#Results in prbi_genepop_2009-11-12_fst.csv and prbi_genepop_2009-11-12_fst.csv.MIG

#Test for Hardy Weinberg
test_HW("PRBI_genepop_2009-11-12.gen.txt", outputFile="prbi_hw_11-12.txt")

#Generate Fst matrix
prbi_11_12_fst <- read.csv("2009-11-12 All Pop All Loci.csv")
prbi_11_12_fst <- as.matrix(prbi_11_12_fst)
prbi_11_12_fst <- prbi_11_12_fst[,-1]

colnames(prbi_11_12_fst) <- c(1, 2, 7, 8, 9, 10, 11, 13, 14, 15, 19, 22)
rownames(prbi_11_12_fst) <- c(1, 2, 7, 8, 9, 10, 11, 13, 14, 15, 19, 22)

##Linearize Fst and Plot Linear Fst vs. Geographic Distance

#Linearize Fst (Fst/(1-Fst))
prbi_11_12_fstlin = as.matrix(prbi_11_12_fst/(1-prbi_11_12_fst))

#Plot Linear Fst vs. Geographic Distance
plot(prbi_geodistance_km, prbi_11_12_fstlin,
     xlab="Pairwise Geographic Distance (km)", ylab="Fst/(1-Fst)", 
     main="2009-11-12 Linearized Fst vs. Geographic Distance", pch=20)

#Add regression to graph
library(smatr)

line_11_12_all <- line.cis(y=c(prbi_11_12_fstlin), x=c(prbi_geodistance_km))
abline(a=line_11_12_all$coef[1], b=line_11_12_all$coef[2], col="red", lwd=1)

#Mantel test

library(vegan)
mantel(prbi_geodistance_km, prbi_11_12_fstlin) 
#r: 0.189 p: 0.11


#2009-11-12 data without population 19:

prbi_no19_11_12_fstlin <- prbi_11_12_fstlin[-11, -11]

plot(prbi_no19_geodist_km, prbi_no19_11_12_fstlin, 
     xlab="Pairwise Geographic Distance (km)", ylab="Fst/(1-Fst)", 
     main="2009-11-12 No Pop 19 Linearized Fst vs. Geographic Distance", pch=20)

#Add regression to graph
line_no19_11_12 <- line.cis(y=c(prbi_no19_11_12_fstlin), x=c(prbi_no19_geodist_km))
abline(a=line_no19_11_12$coef[1], b=line_no19_11_12$coef[2], col="red", lwd=1)

#Mantel test

mantel(prbi_no19_geodist_km, prbi_no19_11_12_fstlin) 
#r: 0.1461 p: 0.207


#2009-11-12 data without populations 13, 14, 15, 22 (East-West transect)

prbi_11_12_EW_fstlin <- prbi_11_12_fstlin[c(-8, -9, -10, -12), c(-8, -9, -10, -12)]

prbi_EW_geodist_km <- prbi_geodistance_km[c(-8, -9, -10, -12), c(-8, -9, -10, -12)]

plot(prbi_EW_geodist_km, prbi_11_12_EW_fstlin, 
     xlab="Pairwise Geographic Distance (km)", ylab="Fst/(1-Fst)", 
     main="2009-11-12 No Pops 13-15, 22 Linearized Fst vs. Geographic Distance", pch=20)

#Add regression to graph
line_11_12_EW <- line.cis(y=c(prbi_11_12_EW_fstlin), x=c(prbi_EW_geodist_km))
abline(a=line_11_12_EW$coef[1], b=line_11_12_EW$coef[2], col="red", lwd=1)

#Mantel test

mantel(prbi_EW_geodist_km, prbi_11_12_EW_fstlin) 
#r: 0.1604 p: 0.235

y <- as.numeric(prbi_11_12_EW_fstlin)
x <- as.numeric(prbi_EW_geodist_km)
mod_EW <- lm(y~x)
summary(mod_EW)
#Multiple R-squared: 0.02574, adjusted R-squared: -0.01174, F-statistic 0.6868 on 1 and 26 DF, p-value: 0.4148
#slope estimate: 1.682e-05, std. error: 2.030e-05, t-value:0.829

#2009-11-12 data without populations 13, 14, 15, 19, 22 (East-West transect wout 19)

prbi_11_12_EW_no19_fstlin <- prbi_11_12_EW_fstlin[-8, -8]

prbi_EW_no19_geodist_km <- prbi_EW_geodist_km[-8, -8]

plot(prbi_EW_no19_geodist_km, prbi_11_12_EW_no19_fstlin, 
     xlab="Pairwise Geographic Distance (km)", ylab="Fst/(1-Fst)", 
     main="2009-11-12 No Pops 13-15, 19, 22 Linearized Fst vs. Geographic Distance", pch=20)

#Add SMA regression to graph
line_11_12_EW_no19 <- line.cis(y=c(prbi_11_12_EW_no19_fstlin), x=c(prbi_EW_no19_geodist_km))
abline(a=line_11_12_EW_no19$coef[1], b=line_11_12_EW_no19$coef[2], col="red", lwd=1)

#Run linear regression
y <- as.numeric(prbi_11_12_EW_no19_fstlin)
x <- as.numeric(prbi_EW_no19_geodist_km)
y
x
class(y)
class(x)

mod_EW_no19 <- lm(y~x)
summary(mod_EW_no19)
#Intercept: Estimate= -2.272e-03, Std. Error= 1.866e-03, t value= -1.218, p value= 0.2383
#x Estimate= 6.253e-05, Std. Error 2.518e-05, t value=2.484, p-value=0.0225

#Multiple R-squared= 0.2451, Adjusted R-squared= 0.2054, F-Statistic=6.168 on 1 and 19 DF,
#p-value=0.0225

coef(mod_EW_no19)
confint(mod_EW_no19)

#Add linear regression to plot
plot(prbi_EW_no19_geodist_km, prbi_11_12_EW_no19_fstlin, 
     xlab="Pairwise Geographic Distance (km)", ylab="Fst/(1-Fst)", 
     main="2009-11-12 No Pops 13-15, 19, 22 Linearized Fst vs. Geographic Distance", pch=20)
abline(mod_EW_no19, col="red")

#Mantel test

mantel(prbi_EW_no19_geodist_km, prbi_11_12_EW_no19_fstlin) 
#r: 0.4951 p: 0.02


#Plot Lin Fst vs. Log Geographic Distance

plot(log(prbi_EW_no19_geodist_km), prbi_11_12_EW_no19_fstlin, 
     xlab="Natural Log of Geographic Distance (km)", ylab="Fst/(1-Fst)", 
     main="2009-11-12 No Pops 13-15, 19, 22 Linearized Fst vs. Log of Geographic Distance", pch=20)

#Run linear regression
y <- as.numeric(prbi_11_12_EW_no19_fstlin)
xlog <- as.numeric(log(prbi_EW_no19_geodist_km))
y
xlog
class(y)
class(xlog)

mod_EW_no19_log <- lm(y~xlog)
summary(mod_EW_no19_log)
#Intercept: Estimate= -0.009870, Std. Error= 0.005673, t value= -1.74, p value= 0.0980
#xlog Estimate= 0.002917, Std. Error 0.001395, t value=2.09, p-value=0.0503

#Multiple R-squared= 0.187, Adjusted R-squared= 0.1442, F-Statistic=4.37 on 1 and 19 DF,
#p-value=0.05026

coef(mod_EW_no19_log)
confint(mod_EW_no19_log)

#Add linear regression to plot
plot(log(prbi_EW_no19_geodist_km), prbi_11_12_EW_no19_fstlin, 
     xlab="Natural Log of Geographic Distance (km)", ylab="Fst/(1-Fst)", 
     main="2009-11-12 No Pops 13-15, 19, 22 Linearized Fst vs. Log of Geographic Distance", pch=20)
abline(mod_EW_no19_log, col="red")

#Mantel test

mantel(log(prbi_EW_no19_geodist_km), prbi_11_12_EW_no19_fstlin) 
#r: 0.4324 p: 0.011



#Calculating Fst with pops 1 and 2 combined as a population

Fst("PRBI_genepop_2009-11-12_comb12.gen.txt", pairs = TRUE, outputFile = "prbi_genepop_2009-11-12_comb12_fst.txt", 
    dataType = "Diploid", verbose = interactive())

#Generate Fst matrix
prbi_11_12_comb12_fst <- read.csv("2009-11-12_comb12_FstMatrix.csv", header=TRUE)
prbi_11_12_comb12_fst <- as.matrix(prbi_11_12_comb12_fst)
prbi_11_12_comb12_fst <- prbi_11_12_comb12_fst[,-1]

colnames(prbi_11_12_comb12_fst) <- c(1, 7, 8, 9, 10, 11, 19)
rownames(prbi_11_12_comb12_fst) <- c(1, 7, 8, 9, 10, 11, 19)

#Linearize Fst (Fst/(1-Fst))
prbi_11_12_comb12_fstlin = as.matrix(prbi_11_12_comb12_fst/(1-prbi_11_12_comb12_fst))





#Mapping the populations
library(maps)
library(mapdata)
#Pull up map of Cebu and Leyte
map(database="world", xlim=c(123, 126), ylim=c(9, 11.5), col="gray", fill=TRUE)
#Add points to the map, representing each population
points(prbi_avglatlong$longitude, prbi_avglatlong$latitude, pch=19, col=1:12)
legend("right", legend=c(1, 2, 7, 8, 9, 10, 11, 13, 14, 15, 19, 22), fill=1:12, col=1:12)

#ggmap requires a google api key
register_google(key="")
library(ggmap)
prbi_map <- get_map(location=c(lon=124, lat=10))



##Trialing different Genepop Fst function inputs

#Trial 1
#pairs=TRUE, the rest are the defaults

Fst("PRBI_genepop_2009-11-12.gen.txt", pairs = TRUE, outputFile = "prbi_genepop_2009-11-12_fst.csv", 
    dataType = "Diploid", verbose = interactive())

#Results in prbi_genepop_2009-11-12_fst.csv and prbi_genepop_2009-11-12_fst.csv.MIG
#Matches "PRBI_11_12_fsts.csv" from GitHub

#Trial 2- all defaults

Fst("PRBI_genepop_2009-11-12.gen.txt",outputFile = "prbi_genepop_defaults_fst.csv")

#With all defaults, it does not generate pairwise matrices

#Trial 3- pairs=TRUE and sizes=TRUE

Fst("PRBI_genepop_2009-11-12.gen.txt", sizes=TRUE, pairs = TRUE, outputFile = "prbi_genepop_sizestrue_fst.csv")
#Does not match "PRBI_11_12_fsts.csv" from GitHub




##Extra/unused:
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





