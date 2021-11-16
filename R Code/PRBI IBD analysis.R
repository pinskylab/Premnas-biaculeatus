#Premnas biaculeatus Genotypes and Analysis

library(genepop)
library(geodist)
library(tidyverse)

##Calculate pairwise Fst between populations, input file from "Downloaded_Data" folder
#pairs=TRUE to get a pairwise Fst matrix

Fst("PRBI_genepop_2009-11-12.gen.txt", pairs = TRUE, outputFile = "prbi_genepop_2009-11-12_fst.csv", 
    dataType = "Diploid", verbose = interactive())

#Results in prbi_genepop_2009-11-12_fst.csv and prbi_genepop_2009-11-12_fst.csv.MIG
#Output file located in "Genepop_Outputs" folder

#Resulting Fst matrix was copied into Microsoft Excel and saved as "2009-11-12 All Pop All Loci.csv" and
#located in "Fst_Matrices" folder


#Test for Hardy Weinberg, input file from "Downloaded_Data" folder
test_HW("PRBI_genepop_2009-11-12.gen.txt", outputFile="prbi_hw_11-12.txt")
#Output file located in "Genepop_Outputs" folder


#Generate Fst matrix
prbi_11_12_fst <- read.csv("2009-11-12 All Pop All Loci.csv")
prbi_11_12_fst <- as.matrix(prbi_11_12_fst)
prbi_11_12_fst <- prbi_11_12_fst[,-1]

colnames(prbi_11_12_fst) <- c(1, 2, 7, 8, 9, 10, 11, 13, 14, 15, 19, 22)
rownames(prbi_11_12_fst) <- c(1, 2, 7, 8, 9, 10, 11, 13, 14, 15, 19, 22)

#To get min and max pairwise Fst
prbi_11_12_fst_short <- prbi_11_12_fst[c(-8, -9, -10, -12), c(-8, -9, -10, -12)]
summary(prbi_11_12_fst_short)

##Linearize Fst and Plot Linear Fst vs. Geographic Distance

#Linearize Fst (Fst/(1-Fst))
prbi_11_12_fstlin = as.matrix(prbi_11_12_fst/(1-prbi_11_12_fst))


##Loading in over water distance matrix
#Calculated using the Google Earth measure tool as the shortest distance between 2 sampling sites
#that does not cross land, distances input into "OverWater_Distance.csv" file located in "Coordinates_Distances" folder

#Load in matrix with over water distances

water_distance <- read.csv("OverWater_Distance.csv")
water_distance <- as.matrix(water_distance)
water_distance <- water_distance[, -1]
colnames(water_distance) <- c(1, 2, 7, 8, 9, 10, 11, 19)
rownames(water_distance) <- c(1, 2, 7, 8, 9, 10, 11, 19)

#Make matrix symmetrical
upperTriangle(water_distance) <- lowerTriangle(water_distance, byrow=TRUE)

library(vegan)

#2009-11-12 data without populations 13, 14, 15, 22 (East-West transect)
#Pops 13, 14, 15, and 22 were excluded due to all having a sample size less than 5 individuals

prbi_11_12_EW_fstlin <- prbi_11_12_fstlin[c(-8, -9, -10, -12), c(-8, -9, -10, -12)]


#Mantel test
mantel(water_distance, prbi_11_12_EW_fstlin)
#r: 0.1469, p: 0.261, 999 permutations 

#Linear regression with over water distance
y <- as.numeric(prbi_11_12_EW_fstlin)
x <- as.numeric(water_distance)
mod_EW_waterdist <- lm(y~x)
summary(mod_EW_waterdist)
#Multiple R-squared: 0.02157, adjusted R-squared: -0.01606, F-statistic 0.5731 on 1 and 26 DF, p: 0.4558
#x estimate: 1.355e-05, std error: 1.789e-05

plot(water_distance, prbi_11_12_EW_fstlin, 
     xlab="Pairwise Geographic Distance (km)", ylab="Fst/(1-Fst)")
abline(mod_EW_waterdist, col="red")


#Excluding sites 13, 14, 15, 19, and 22, as 19 does not appear to follow an IBD pattern when plotted

prbi_11_12_EW_no19_fstlin <- prbi_11_12_EW_fstlin[-8, -8]

water_distance_no19 <- water_distance[-8, -8]

#Mantel test with over water distance
mantel(water_distance_no19, prbi_11_12_EW_no19_fstlin)
#r: 0.4563, p: 0.03, 5039 permutations 

#Linear regression with over water distance
y <- as.numeric(prbi_11_12_EW_no19_fstlin)
x <- as.numeric(water_distance_no19)

mod_EW_no19_waterdist <- lm(y~x)
summary(mod_EW_no19_waterdist)
#Intercept estimate: -1.791e-03, std. error: 1.274e-03, t-value:-1.406, p:0.16751
#x estimate: 5.393e-05, std error: 1.663e-05, t-value: 3.243, p:0.00239
#Multiple R-squared: 0.2082, adjusted R-squared: 0.1884, F-statistic 10.52 on 1 and 40 DF, p: 0.002388

plot(water_distance_no19, prbi_11_12_EW_no19_fstlin, 
     xlab="Pairwise Geographic Distance (km)", ylab="Fst/(1-Fst)")
abline(mod_EW_no19_waterdist, col="red")


#Plotting geographic distance vs. genetic distance, showing pop 19 in circles

distanceframe <- tibble(GeneticDistance=as.vector(prbi_11_12_EW_no19_fstlin), 
                        WaterDistance=as.vector(water_distance_no19))
distanceframe <- drop_na(distanceframe)

pop19 <- tibble(WaterDistance = as.vector(water_distance[8,]),
                GeneticDistance = as.vector(prbi_11_12_EW_fstlin[8,]))
pop19 <- drop_na(pop19)

ggplot(data=distanceframe, aes(x=distanceframe$WaterDistance, 
                               y=distanceframe$GeneticDistance)) +
    geom_point(size=2.5) +
    geom_point(data=pop19, aes(x=WaterDistance, y=GeneticDistance), size=3, shape=1) +
    geom_smooth(method="lm", color="black") +
    theme_bw() +
    xlab("Geographic Distance (km)") +
    ylab("Genetic Distance (Fst/(1-Fst)")



#Calculating Fst with pops 1 and 2 combined as a population
#This was needed for the potential connectivity analyses, as sites 1 and 2 were too close
#together for two distinct particle release sites to be identified
#Input file from "Data" folder

Fst("PRBI_genepop_2009-11-12_comb12.gen.txt", pairs = TRUE, outputFile = "prbi_genepop_2009-11-12_comb12_fst.txt", 
    dataType = "Diploid", verbose = interactive())
#Output file located in "Genepop_Outputs" folder

#Resulting matrix copied into Microsoft Excel and saved as "2009-11-12_comb12_FstMatrix.csv" and located in
#"Fst_Matrices" folder

#Generate Fst matrix
prbi_11_12_comb12_fst <- read.csv("2009-11-12_comb12_FstMatrix.csv", header=TRUE)
prbi_11_12_comb12_fst <- as.matrix(prbi_11_12_comb12_fst)
prbi_11_12_comb12_fst <- prbi_11_12_comb12_fst[,-1]

colnames(prbi_11_12_comb12_fst) <- c(1, 7, 8, 9, 10, 11, 19)
rownames(prbi_11_12_comb12_fst) <- c(1, 7, 8, 9, 10, 11, 19)

#Linearize Fst (Fst/(1-Fst))
prbi_11_12_comb12_fstlin = as.matrix(prbi_11_12_comb12_fst/(1-prbi_11_12_comb12_fst))






