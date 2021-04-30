##Investigating how Ocean Currents may Impact PRBI Dispersal##

library(vegan) #Mantel test
library(gdata) #Has functions for transposing a symmetrical matrix

#Load in connectivity matrix

ocean_matrix <- read.csv("Oceanographic_Matrix.csv")
ocean_matrix <- as.matrix(ocean_matrix)
ocean_matrix <- ocean_matrix[, -1]
colnames(ocean_matrix) <- c(1, 2, 7, 8, 9, 10, 11, 19)
rownames(ocean_matrix) <- c(1, 2, 7, 8, 9, 10, 11, 19)

#Mantel test of connectivity matrix (oceanographic distance) vs. genetic distance (linearized Fst)

#Load in Fst matrix for all populations
prbi_11_12_fst_noNA <- read.csv("2009-11-12 All Pop No NA Fst Matrix.csv")

prbi_11_12_fst_noNA <- as.matrix(prbi_11_12_fst_noNA)
prbi_11_12_fst_noNA <- prbi_11_12_fst_noNA[,-1]

colnames(prbi_11_12_fst_noNA) <- c(1, 2, 7, 8, 9, 10, 11, 13, 14, 15, 19, 22)
rownames(prbi_11_12_fst_noNA) <- c(1, 2, 7, 8, 9, 10, 11, 13, 14, 15, 19, 22)

##Linearize Fst

#Linearize Fst (Fst/(1-Fst))
prbi_11_12_noNA_fstlin = as.matrix(prbi_11_12_fst_noNA/(1-prbi_11_12_fst_noNA), na.rm=TRUE)
#Linear Fst matrix, has lower triangle filled in but otherwise is NAs

#Fill in the values on the other side of the diagonal (this is a symmetric matrix, so the values are mirrored across the diagonal)

upperTriangle(prbi_11_12_noNA_fstlin) <- lowerTriangle(prbi_11_12_noNA_fstlin, byrow=TRUE) #This works to get symmetrical matrix

prbi_11_12_noNA_fstlin

#Now delete columns and rows corresponding to Pops 13, 14, 15, and 22- they are being excluded due to sample size lower than 10

prbi_11_12_EW_full_fstlin <- prbi_11_12_noNA_fstlin[c(-8, -9, -10, -12), c(-8, -9, -10, -12)]
prbi_11_12_EW_full_fstlin

#Perform Mantel test to assess correlation between connectivity matrix (gives probability of dispersal between populations) and
#genetic distance (linearized Fst)

mantel(ocean_matrix, prbi_11_12_EW_full_fstlin)
# r: -0.4342, p:0.989
#Mantel statistic based on Pearson's product-moment correlation, 999 permutations, free permutation

#Run linear regression
y <- as.numeric(prbi_11_12_EW_full_fstlin)
x <- as.numeric(ocean_matrix)
y
x
class(y)
class(x)

mod_EW_ocean <- lm(y~x)
summary(mod_EW_ocean)
#adjusted R-squared:0.02729 p: 0.1166
#x estimate: -0.0390143, std. error: 0.0244654

#Plot data
plot(ocean_matrix, prbi_11_12_EW_full_fstlin, 
     xlab="Probability of Larval Dispersal between Populations", ylab="Fst/(1-Fst)", 
      pch=20)
abline(mod_EW_ocean, col="red")

#Run Mantel test with lower triangle of oceanographic matrix
lower_ocean_matrix <- read.csv("Oceanographic_Matrix.csv")
lower_ocean_matrix <- as.matrix(lower_ocean_matrix)
lower_ocean_matrix <- lower_ocean_matrix[, -1]
colnames(lower_ocean_matrix) <- c(1, 2, 7, 8, 9, 10, 11, 19)
rownames(lower_ocean_matrix) <- c(1, 2, 7, 8, 9, 10, 11, 19)
lower_ocean_matrix[upper.tri(lower_ocean_matrix, diag=T)] = NA

mantel(lower_ocean_matrix, prbi_11_12_EW_fstlin)
# r: -0.4342, p:0.984
#Mantel statistic based on Pearson's product-moment correlation, 999 permutations, free permutation

#Run Mantel test with upper triangle of oceanographic matrix
#This code looks a little counter-intuitive because the mantel function takes the values from the lower triangle of a square
#matrix, so I need to move the dispersal probabilities from the upper triangle of the ocean_matrix to the lower triangle

upper_ocean_matrix <- as.matrix(ocean_matrix)
upper_ocean_matrix[lower.tri(upper_ocean_matrix, diag=T)] = NA
lowerTriangle(upper_ocean_matrix) <- upperTriangle(upper_ocean_matrix, byrow=TRUE)
upper_ocean_matrix[upper.tri(upper_ocean_matrix, diag=T)] = NA

mantel(upper_ocean_matrix, prbi_11_12_EW_fstlin)
# r:0.2851, p:0.085
#Mantel statistic based on Pearson's product-moment correlation, 999 permutations, free permutation

##Run partial Mantel test of oceanographic distance vs. genetic distance, while controlling for geographic distance

#Load in matrix with over water distances

water_distance <- read.csv("OverWater_Distance.csv")
water_distance <- as.matrix(water_distance)
water_distance <- water_distance[, -1]
colnames(water_distance) <- c(1, 2, 7, 8, 9, 10, 11, 19)
rownames(water_distance) <- c(1, 2, 7, 8, 9, 10, 11, 19)

#Make matrix symmetrical
upperTriangle(water_distance) <- lowerTriangle(water_distance, byrow=TRUE)

#Partial Mantel test controlling for over water geographic distance
#This is based on the lower triangle of the ocean matrix
mantel.partial(ocean_matrix, prbi_11_12_EW_full_fstlin, water_distance)
#Partial Mantel statistic r: -0.4157, p: 0.978

#Partial Mantel using the upper triangle of the ocean matrix

mantel.partial(upper_ocean_matrix, prbi_11_12_EW_fstlin, water_distance)
#r: 0.3004, p: 0.082

