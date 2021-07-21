## Potential Connectivity Matrices
#From Diane Thompson- for 10-day PLD

##Loading in matrices with tidync package

install.packages("tidync")
library(tidync)

library(ncdf4)
connect_matrix <- nc_open("CORAL-connect_25_historical_pconnect_10-day_noedges2.nc")
print(connect_matrix)  #This gives the meta-data for the matrix
#[dst_site, src_site] - dst=destination site, src=source site

connect_data <- ncvar_get(connect_matrix, varid="pconnect") #this is the actual potential connectivity matrix
print(connect_data)
dim(connect_data)
head(connect_data)

latlon_matrix <- nc_open("coral_all_release_site_coarse.nc")
print(latlon_matrix)

lon_data <- ncvar_get(latlon_matrix, varid="lon")
print(lon_data)
dim(lon_data)
lon_data[1:10, 1:10]

lat_data <- ncvar_get(latlon_matrix, varid="lat")
print(lat_data)
lat_data[1:10, 1:10]

release_site_data <- ncvar_get(latlon_matrix, varid="release_site")
print(release_site_data)
dim(release_site_data)
release_site_data[1:10, 1:10]
sort(unique(c(release_site_data)))
range(unique(c(release_site_data))) #2947
i <- which(c(release_site_data)==1)  #setting i as release site=1,
#lumped coordinates into 2947 release sites

c(lat_data)[i]  #lat of release site 1
c(lon_data)[i]  #lon of release site 1

plot(c(lon_data)[i], c(lat_data)[i]) 

which.min(abs(c(lat_data) - MYSITELAT)) #use this to find a release site that matches
#each coordinate

#Possible release sites that match my coordinates: 1925-1929; 1955-1959; 1993, 1994, 1997; 2041-2042
#1925-1929, 1955-1959, 1993, 1994, 1997, 2041-2042 do not match

i <- which(c(lat_data) > 9) 
c(lat_data)[i]  
c(lon_data)[i]
pos_coord <- c(release_site_data)[i]
View(pos_coord) #vector with release sites where latitude is greater than 9

i <- which(c(lon_data) > 123)
pos_lon <- c(release_site_data)[i]
View(pos_lon)
pos_lon <- fun.zero.omit(pos_lon)
pos_lon <- sort(pos_lon)
pos_lon <- table(pos_lon)

library(GLDEX)
pos_coord <- fun.zero.omit(pos_coord) #omits zeros from the vector
pos_coord <- sort(pos_coord) #orders from least to greatest
pos_coord <- table(pos_coord) #gives frequency table, so each release site is only listed once

i <- which(c(release_site_data)==2481)
c(lat_data)[i]  
c(lon_data)[i]

#Possible matching release sites: 2375, 2401 (pop 7?), 2402, 2426, 2427, 2428, 2454, 2455, 2457, 
#2458, 2479, 2480, 2481, 2483 (pop 19?), 2484

#Pop 1- site 2480  (site 2481 is also similar)
#Pop 2- site 2480
#Pop 7- site 2401
#Pop 8- site 2426
#Pop 9- site 2427 (maybe- is a little south)
#Pop 10- site 2454 (maybe)
#Pop 11- site 2455 (maybe)
#Pop 19- site 2483

#Make vector with pconnect between pop 19 and other pops

#1st in direction of pop 19 to other pops
#[dst site, src site] would be [other pop, pop 19], so [other site, 2483]
#Order of vector: Pop 1, 2, 7, 8, 9, 10, 11, 19 
pconnect_pop19src <- c(connect_data[2480, 2483], connect_data[2480, 2483], connect_data[2401, 2483], connect_data[2426, 2483],
                      connect_data[2427, 2483], connect_data[2454, 2483], connect_data[2455, 2483])
pconnect_pop19src

connect_data[2483, 2483] #probability of self-recruitment?

#2nd in direction of other pops to pop 19
pconnect_pop19dst <- c(connect_data[2483, 2480], connect_data[2483, 2480], connect_data[2483, 2401], connect_data[2483, 2426],
                       connect_data[2483, 2427], connect_data[2483, 2454], connect_data[2483, 2455])

pconnect_pop19dst

pconnect_avg <- apply(rbind(pconnect_pop19src, pconnect_pop19dst), 2, mean)
pconnect_avg

pop_numbers <- c(1, 2, 7, 8, 9, 10, 11)

pconnect_pop19src_df <- data.frame(pop_numbers, pconnect_pop19src)

pconnect_pop19dst_df <- data.frame(pop_numbers, pconnect_pop19dst)

pconnect_avg_df <- data.frame(pop_numbers, pconnect_avg)

#probability of pops 7, 8, and 9 dispersing passively to pop 19 is 0, whereas probability of pop 19
#dispersing to pops 7, 8, and 9 is higher than others.

#Then plot Fst on y-axis and probability of connection on the x-axis, for now average the two
#potential connectivity values that come up (the matrix isn't symmetrical so there will be 2 values
#between each population)

pop19 #Data frame with pairwise linear Fst between pop 19 and the other pops
pop19_dist <- c(0,0)
pop19 <- rbind(pop19, pop19_dist)
pop19 #Adding that the geographic and genetic distance from pop 19 to pop 19 is 0
pop19 <- pop19[-c(8),] #deleting row with pop 19

pconnect_avg_df <- data.frame(pop_numbers, pconnect_avg, pop19$GeneticDistance)
pconnect_avg_df

#Plotting average pconnect vs. lin Fst
mod_pconnect_avg <- lm(pconnect_avg_df$pop19.GeneticDistance ~ pconnect_avg_df$pconnect_avg)
summary(mod_pconnect_avg) #Adjusted R-squared:0.823, p:0.003

plot(pconnect_avg_df$pconnect_avg, pop19$GeneticDistance, pch=19)
abline(mod_pconnect_avg, col="red")

pconnect_avg_df$pop_numbers <- as.factor(pconnect_avg_df$pop_numbers)
class(pconnect_avg_df$pop_numbers)

ggplot(data=pconnect_avg_df, aes(x=pconnect_avg, 
                                             y=pop19.GeneticDistance, color=pop_numbers)) +
  geom_point(size=1.75) +
  geom_smooth(method="lm", se=FALSE) +
  theme_bw() +
  xlab("PConnect") +
  ylab("Fst/(1-Fst)") +
  scale_color_manual(values = c("1" = "purple",
                                "2"="orange",
                                "7"="steelblue",
                                "8"="red",
                                "9"="gold",
                                "10"="green",
                                "11"="pink")) 

ggplot(data=pconnect_avg_df, aes(x=pconnect_avg, 
                                 y=pop19.GeneticDistance, color=pop_numbers)) +
  geom_point(size=2.5) +
  geom_smooth(method="lm", se=FALSE, color="steelblue") +
  theme_bw() +
  xlab("PConnect") +
  ylab("Fst/(1-Fst)")

pconnect_avg_df


#Plotting pconnect with pop 19 as the source vs. lin Fst

mod_pconnect_src <- lm(pop19$GeneticDistance ~ pconnect_pop19src_df$pconnect_pop19src)
summary(mod_pconnect_src) #Adjusted R-squared: 0.8255, p:0.002895

plot(pconnect_pop19src_df$pconnect_pop19src, pop19$GeneticDistance, pch=19)
abline(mod_pconnect_src, col="red")


#Plotting pconnect with pop 19 as the destination vs. lin Fst

mod_pconnect_dst <- lm(pop19$GeneticDistance ~ pconnect_pop19dst_df$pconnect_pop19dst)
summary(mod_pconnect_dst) #Adjusted R-squared: -0.1227, p:0.5828

plot(pconnect_pop19dst_df$pconnect_pop19dst, pop19$GeneticDistance, pch=19)
abline(mod_pconnect_dst, col="red")  #pop 19 as the destination looks very different,
#has pop 11 with a higher probability of connection than population 1 or 2- look into this


pconnect_pop19src_df <- data.frame(pop_numbers, pconnect_pop19src, pop19$GeneticDistance)

pconnect_pop19dst_df <- data.frame(pop_numbers, pconnect_pop19dst, pop19$GeneticDistance)

pconnect_pop19src_df$pop_numbers <- as.factor(pconnect_pop19src_df$pop_numbers)
pconnect_pop19dst_df$pop_numbers <- as.factor(pconnect_pop19dst_df$pop_numbers)

#Ggplot with pop 19 as the source
ggplot(data=pconnect_pop19src_df, aes(x=pconnect_pop19src, 
                                 y=pop19.GeneticDistance, color=pop_numbers)) +
  geom_point(size=2.5) +
  geom_smooth(method="lm", se=FALSE, color="steelblue") +
  theme_bw() +
  xlab("PConnect- Pop 19 Source") +
  ylab("Fst/(1-Fst)")

#Ggplot with pop 19 as the destination
ggplot(data=pconnect_pop19dst_df, aes(x=pconnect_pop19dst, 
                                      y=pop19.GeneticDistance, color=pop_numbers)) +
  geom_point(size=2.5) +
  theme_bw() +
  xlab("PConnect- Pop 19 Destination") +
  ylab("Fst/(1-Fst)")

#Could map coordinates used from the pconnect matrix and compare to our populations, 
#add arrows to map weighted based on pconnect


#Create pairwise pconnect matrix

#Pop 1- site 2480  
#Pop 2- site 2480
#Pop 7- site 2401
#Pop 8- site 2426
#Pop 9- site 2427 
#Pop 10- site 2454 
#Pop 11- site 2455 
#Pop 19- site 2483

#Lower triangle pconnect matrix 
#Pop 19 as the destination

pconnect_lowertri_vector <- c(NA, NA, NA, NA, NA, NA, NA, NA,
                              connect_data[2480, 2480], NA, NA, NA, NA, NA, NA, NA,
                              connect_data[2401, 2480], connect_data[2401, 2480], NA, NA, NA, NA, NA, NA,
                              connect_data[2426, 2480], connect_data[2426, 2480], connect_data[2426, 2401], NA, NA, NA, NA, NA,
                              connect_data[2427, 2480], connect_data[2427, 2480], connect_data[2427, 2401], connect_data[2427, 2426], NA, NA, NA, NA,
                              connect_data[2454, 2480], connect_data[2454, 2480], connect_data[2454, 2401], connect_data[2454, 2426], connect_data[2454, 2427], NA, NA, NA,
                              connect_data[2455, 2480], connect_data[2455, 2480], connect_data[2455, 2401], connect_data[2455, 2426], connect_data[2455, 2427], connect_data[2455, 2454], NA, NA,
                              connect_data[2483, 2480], connect_data[2483, 2480], connect_data[2483, 2401], connect_data[2483, 2426], connect_data[2483, 2427], connect_data[2483, 2454], connect_data[2483, 2455], NA)
pconnect_lowertri <- matrix(pconnect_lowertri_vector, nrow=8)
pconnect_lowertri
colnames(pconnect_lowertri) <- c(1, 2, 7, 8, 9, 10, 11, 19)
rownames(pconnect_lowertri) <- c(1, 2, 7, 8, 9, 10, 11, 19)

#Upper triangle pconnect matrix
#Pop 19 as the source

pconnect_uppertri_vector <- c(NA, NA, NA, NA, NA, NA, NA, NA,
                              connect_data[2480, 2480], NA, NA, NA, NA, NA, NA, NA,
                              connect_data[2480, 2401], connect_data[2480, 2401], NA, NA, NA, NA, NA, NA,
                              connect_data[2480, 2426], connect_data[2480, 2426], connect_data[2401, 2426], NA, NA, NA, NA, NA,
                              connect_data[2480, 2427], connect_data[2480, 2427], connect_data[2401, 2427], connect_data[2426, 2427], NA, NA, NA, NA,
                              connect_data[2480, 2454], connect_data[2480, 2454], connect_data[2401, 2454], connect_data[2426, 2454], connect_data[2427, 2454], NA, NA, NA,
                              connect_data[2480, 2455], connect_data[2480, 2455], connect_data[2401, 2455], connect_data[2426, 2455], connect_data[2427, 2455], connect_data[2454, 2455], NA, NA,
                           connect_data[2480, 2483], connect_data[2480, 2483], connect_data[2401, 2483], connect_data[2426, 2483], connect_data[2427, 2483], connect_data[2454, 2483], connect_data[2455, 2483], NA)

pconnect_uppertri <- matrix(pconnect_uppertri_vector, nrow=8)
pconnect_uppertri
colnames(pconnect_uppertri) <- c(1, 2, 7, 8, 9, 10, 11, 19)
rownames(pconnect_uppertri) <- c(1, 2, 7, 8, 9, 10, 11, 19)

#Compute average pconnect matrix

pconnect_avg_matrix <- apply(rbind(pconnect_lowertri_vector, pconnect_uppertri_vector), 2, mean)
pconnect_avg_matrix <- matrix(pconnect_avg_matrix, nrow=8)
pconnect_avg_matrix
colnames(pconnect_avg_matrix) <- c(1, 2, 7, 8, 9, 10, 11, 19)
rownames(pconnect_avg_matrix) <- c(1, 2, 7, 8, 9, 10, 11, 19)

#Flip values across diagonal, so they are in the lower triangle of the matrix

library(gdata)
lowerTriangle(pconnect_avg_matrix) <- upperTriangle(pconnect_avg_matrix, byrow=TRUE)

#Mantel test between average pconnect and linearized Fst 

library(vegan)
mantel(pconnect_avg_matrix, prbi_11_12_EW_full_fstlin) 
#Mantel statistic r: -0.2447, p-value: 0.913, 999 permutations

plot(pconnect_avg_matrix, prbi_11_12_EW_full_fstlin, 
     xlab="Probability of Larval Dispersal between Populations", ylab="Fst/(1-Fst)", 
     pch=20)

#Run a linear regression between average pconnect and linearized Fst
y <- as.numeric(prbi_11_12_EW_full_fstlin)
x <- as.numeric(pconnect_avg_matrix)

mod_pconnect_matrix <- lm(y~x)
summary(mod_pconnect_matrix)
#Adjusted R-squared: 0.04248, p-value: 0.06909


plot(pconnect_avg_matrix, prbi_11_12_EW_full_fstlin, 
     xlab="Probability of Larval Dispersal between Populations", ylab="Fst/(1-Fst)", 
     pch=20)
abline(mod_pconnect_matrix, col="red")

#Use performance package to see how outliers are influencing model
library(performance)
check_model(mod_pconnect_matrix)
plot(mod_pconnect_matrix)

#Looks like the 0.3 pconnect point has too much leverage, this is the pconnect between 
#pops 1 and 2 which are the same release site, at this point the pconnect of self-recruitment is not
#included, so I'll try excluding it, may need to run again while including self-recruitment
#0.105 pconnect is between pops 9 and 8
#0.127 pconnect is between pops 8 and 10- geographic distance seems to have the better 
#correlation her as pops 8 and 9 have a lin Fst closer to 0 than pops 8 and 10 despite the 
#higher pconnect
#0.0777 pconnect is between pops 9 and 10

#Remove pconnect between sites 1 and 2 (self-recruitment pconnect)

pconnect_avg_matrix_noself <- pconnect_avg_matrix
pconnect_avg_matrix_noself[1, 2] <- NA
pconnect_avg_matrix_noself[2,1] <- NA

prbi_11_12_EW_full_fstlin_noself <- prbi_11_12_EW_full_fstlin
prbi_11_12_EW_full_fstlin_noself[1, 2] <- NA
prbi_11_12_EW_full_fstlin_noself[2,1] <- NA

mantel(pconnect_avg_matrix_noself, prbi_11_12_EW_full_fstlin_noself)
#Won't run because missing observations

y2 <- as.numeric(prbi_11_12_EW_full_fstlin_noself)
x2 <- as.numeric(pconnect_avg_matrix_noself)

mod_pconnect_matrix_noself <- lm(y2~x2)
summary(mod_pconnect_matrix_noself)
#Adjusted R-squared: 0.2216, p-value: 0.000194


plot(pconnect_avg_matrix_noself, prbi_11_12_EW_full_fstlin_noself, 
     xlab="Probability of Larval Dispersal between Populations", ylab="Fst/(1-Fst)", 
     pch=20)
abline(mod_pconnect_matrix_noself, col="red")

check_model(mod_pconnect_matrix_noself)
plot(mod_pconnect_matrix_noself)
#Residuals vs. leverage plot looks better than before, but still not great


##Partial Mantel test of oceanographic distance vs. genetic distance, while controlling for geographic distance (in km)

#Using over-water geographic distances calculated using ruler tool in Google Earth

water_distance <- read.csv("OverWater_Distance.csv")
water_distance <- as.matrix(water_distance)
water_distance <- water_distance[, -1]
colnames(water_distance) <- c(1, 2, 7, 8, 9, 10, 11, 19)
rownames(water_distance) <- c(1, 2, 7, 8, 9, 10, 11, 19)

#Make matrix symmetrical
upperTriangle(water_distance) <- lowerTriangle(water_distance, byrow=TRUE)

#Partial Mantel test controlling for over water geographic distance

mantel.partial(pconnect_avg_matrix, prbi_11_12_EW_full_fstlin, water_distance)
#Mantel statistic r: -0.2004, p-value: 0.846, 999 permutations


##Changing Pop 1 to release site 2481##

#Pop 1- site 2481
#Pop 2- site 2480
#Pop 7- site 2401
#Pop 8- site 2426
#Pop 9- site 2427 
#Pop 10- site 2454
#Pop 11- site 2455
#Pop 19- site 2483

#Lower triangle pconnect matrix 
#Pop 19 as the destination

pconnect_lowertri_2481_vector <- c(NA, NA, NA, NA, NA, NA, NA, NA,
                              connect_data[2480, 2481], NA, NA, NA, NA, NA, NA, NA,
                              connect_data[2401, 2481], connect_data[2401, 2480], NA, NA, NA, NA, NA, NA,
                              connect_data[2426, 2481], connect_data[2426, 2480], connect_data[2426, 2401], NA, NA, NA, NA, NA,
                              connect_data[2427, 2481], connect_data[2427, 2480], connect_data[2427, 2401], connect_data[2427, 2426], NA, NA, NA, NA,
                              connect_data[2454, 2481], connect_data[2454, 2480], connect_data[2454, 2401], connect_data[2454, 2426], connect_data[2454, 2427], NA, NA, NA,
                              connect_data[2455, 2481], connect_data[2455, 2480], connect_data[2455, 2401], connect_data[2455, 2426], connect_data[2455, 2427], connect_data[2455, 2454], NA, NA,
                              connect_data[2483, 2481], connect_data[2483, 2480], connect_data[2483, 2401], connect_data[2483, 2426], connect_data[2483, 2427], connect_data[2483, 2454], connect_data[2483, 2455], NA)
pconnect_lowertri_2481 <- matrix(pconnect_lowertri_2481_vector, nrow=8)
pconnect_lowertri_2481
colnames(pconnect_lowertri_2481) <- c(1, 2, 7, 8, 9, 10, 11, 19)
rownames(pconnect_lowertri_2481) <- c(1, 2, 7, 8, 9, 10, 11, 19)

#Upper triangle pconnect matrix
#Pop 19 as the source

pconnect_uppertri_2481_vector <- c(NA, NA, NA, NA, NA, NA, NA, NA,
                              connect_data[2481, 2480], NA, NA, NA, NA, NA, NA, NA,
                              connect_data[2481, 2401], connect_data[2480, 2401], NA, NA, NA, NA, NA, NA,
                              connect_data[2481, 2426], connect_data[2480, 2426], connect_data[2401, 2426], NA, NA, NA, NA, NA,
                              connect_data[2481, 2427], connect_data[2480, 2427], connect_data[2401, 2427], connect_data[2426, 2427], NA, NA, NA, NA,
                              connect_data[2481, 2454], connect_data[2480, 2454], connect_data[2401, 2454], connect_data[2426, 2454], connect_data[2427, 2454], NA, NA, NA,
                              connect_data[2481, 2455], connect_data[2480, 2455], connect_data[2401, 2455], connect_data[2426, 2455], connect_data[2427, 2455], connect_data[2454, 2455], NA, NA,
                              connect_data[2481, 2483], connect_data[2480, 2483], connect_data[2401, 2483], connect_data[2426, 2483], connect_data[2427, 2483], connect_data[2454, 2483], connect_data[2455, 2483], NA)

pconnect_uppertri_2481 <- matrix(pconnect_uppertri_2481_vector, nrow=8)
pconnect_uppertri_2481
colnames(pconnect_uppertri_2481) <- c(1, 2, 7, 8, 9, 10, 11, 19)
rownames(pconnect_uppertri_2481) <- c(1, 2, 7, 8, 9, 10, 11, 19)

#Compute average pconnect matrix

pconnect_avg_2481_matrix <- apply(rbind(pconnect_lowertri_2481_vector, pconnect_uppertri_2481_vector), 2, mean)
pconnect_avg_2481_matrix <- matrix(pconnect_avg_2481_matrix, nrow=8)
pconnect_avg_2481_matrix
colnames(pconnect_avg_2481_matrix) <- c(1, 2, 7, 8, 9, 10, 11, 19)
rownames(pconnect_avg_2481_matrix) <- c(1, 2, 7, 8, 9, 10, 11, 19)

#Flip values across diagonal, so they are in the lower triangle of the matrix

library(gdata)
lowerTriangle(pconnect_avg_2481_matrix) <- upperTriangle(pconnect_avg_2481_matrix, byrow=TRUE)

#Mantel test between average pconnect and linearized Fst 

library(vegan)
mantel(pconnect_avg_2481_matrix, prbi_11_12_EW_full_fstlin) 
#Mantel statistic r: -0.4809, p-value: 0.994, 999 permutations

plot(pconnect_avg_2481_matrix, prbi_11_12_EW_full_fstlin, 
     xlab="Probability of Larval Dispersal between Populations", ylab="Fst/(1-Fst)", 
     pch=20)

#Run a linear regression between average pconnect and linearized Fst
y <- as.numeric(prbi_11_12_EW_full_fstlin)
x <- as.numeric(pconnect_avg_2481_matrix)

mod_pconnect_2481_matrix <- lm(y~x)
summary(mod_pconnect_2481_matrix)
#Adjusted R-squared: 0.217, p-value: 0.0001758

cor.test(x, y, method='pearson')
#cor: -0.480878, p-value: 0.0001758


plot(pconnect_avg_2481_matrix, prbi_11_12_EW_full_fstlin, 
     xlab="Probability of Larval Dispersal between Populations", ylab="Fst/(1-Fst)", 
     pch=20)
abline(mod_pconnect_2481_matrix, col="red")

#Partial Mantel test controlling for over water geographic distance

mantel.partial(pconnect_avg_2481_matrix, prbi_11_12_EW_full_fstlin, water_distance)
#Mantel statistic r: -0.4729, p-value: 0.996, 999 permutations


##Plot Pop 19 compared to other pops with 2481 as pop 1 release site##

#Make vector with pconnect between pop 19 and other pops

#1st in direction of pop 19 to other pops
#[dst site, src site] would be [other pop, pop 19], so [other site, 2483]
#Order of vector: Pop 1, 2, 7, 8, 9, 10, 11, 19 
pconnect_pop19src_2481 <- c(connect_data[2481, 2483], connect_data[2480, 2483], connect_data[2401, 2483], connect_data[2426, 2483],
                       connect_data[2427, 2483], connect_data[2454, 2483], connect_data[2455, 2483])
pconnect_pop19src_2481

#2nd in direction of other pops to pop 19
pconnect_pop19dst_2481 <- c(connect_data[2483, 2481], connect_data[2483, 2480], connect_data[2483, 2401], connect_data[2483, 2426],
                       connect_data[2483, 2427], connect_data[2483, 2454], connect_data[2483, 2455])

pconnect_pop19dst_2481

pconnect_avg_2481 <- apply(rbind(pconnect_pop19src_2481, pconnect_pop19dst_2481), 2, mean)
pconnect_avg_2481

pop_numbers <- c(1, 2, 7, 8, 9, 10, 11)

pconnect_avg_2481_df <- data.frame(pop_numbers, pconnect_avg_2481, pop19$GeneticDistance)
pconnect_avg_2481_df

#Plotting average pconnect vs. lin Fst
mod_pconnect_2481_avg <- lm(pconnect_avg_2481_df$pop19.GeneticDistance ~ pconnect_avg_2481_df$pconnect_avg_2481)
summary(mod_pconnect_2481_avg) #Adjusted R-squared:0.7537, p:0.00702

cor.test(pconnect_avg_2481_df$pconnect_avg_2481, pconnect_avg_2481_df$pop19.GeneticDistance, method = 'pearson')
#cor: -0.89, p-value: 0.007, df=5

plot(pconnect_avg_2481_df$pconnect_avg_2481, pop19$GeneticDistance, pch=19)
abline(mod_pconnect_2481_avg, col="red")

pconnect_avg_2481_df$pop_numbers <- as.factor(pconnect_avg_2481_df$pop_numbers)
class(pconnect_avg_df$pop_numbers)

ggplot(data=pconnect_avg_2481_df, aes(x=pconnect_avg_2481, 
                                 y=pop19.GeneticDistance, color=pop_numbers)) +
  geom_point(size=2.5) +
  geom_smooth(method="lm", se=FALSE, color="steelblue") +
  theme_bw() +
  xlab("PConnect") +
  ylab("Fst/(1-Fst)")


#Excluding pop 19 from pairwise matrices while using release site 2481 for pop 1#

pconnect_avg_2481_no19_matrix <- pconnect_avg_2481_matrix[-8, -8]

upperTriangle(prbi_11_12_EW_no19_fstlin) <- lowerTriangle(prbi_11_12_EW_no19_fstlin, byrow=TRUE)

#Mantel test between average pconnect and linearized Fst 

mantel(pconnect_avg_2481_no19_matrix, prbi_11_12_EW_no19_fstlin) 
#Mantel statistic r: -0.4296, p-value: 0.967, 5039 permutations

plot(pconnect_avg_2481_no19_matrix, prbi_11_12_EW_no19_fstlin, 
     xlab="Probability of Larval Dispersal between Populations", ylab="Fst/(1-Fst)", 
     pch=20)

#Run a linear regression between average pconnect and linearized Fst
y <- as.numeric(prbi_11_12_EW_no19_fstlin)
x <- as.numeric(pconnect_avg_2481_no19_matrix)

mod_pconnect_2481_no19_matrix <- lm(y~x)
summary(mod_pconnect_2481_no19_matrix)
#Adjusted R-squared: 0.1642, p-value: 0.004524

cor.test(x, y, method='pearson')
#cor: -0.4295782, p-value: 0.004524


plot(pconnect_avg_2481_no19_matrix, prbi_11_12_EW_no19_fstlin, 
     xlab="Probability of Larval Dispersal between Populations", ylab="Fst/(1-Fst)", 
     pch=20)
abline(mod_pconnect_2481_no19_matrix, col="red")

#Partial Mantel test controlling for over water geographic distance

water_distance_no19 <- water_distance[-8, -8]

mantel.partial(pconnect_avg_2481_no19_matrix, prbi_11_12_EW_no19_fstlin, water_distance_no19)
#Mantel statistic r: -0.2141, p-value: 0.776, 5039 permutations


##Pairwise matrix with pops 1 and 2 Fst combined and release site 2480 for pop 1#
#Pop 1- site 2480
#Pop 7- site 2401
#Pop 8- site 2426
#Pop 9- site 2427
#Pop 10- site 2454 
#Pop 11- site 2455 
#Pop 19- site 2483

pconnect_avg_comb12_matrix <- pconnect_avg_matrix[-2, -2]

upperTriangle(prbi_11_12_comb12_fstlin) <- lowerTriangle(prbi_11_12_comb12_fstlin, byrow=TRUE)

#Mantel test between average pconnect and linearized Fst 

mantel(pconnect_avg_comb12_matrix, prbi_11_12_comb12_fstlin) 
#Mantel statistic r: -0.4417, p-value: 0.975, 5039 permutations

plot(pconnect_avg_comb12_matrix, prbi_11_12_comb12_fstlin, 
     xlab="Probability of Larval Dispersal between Populations", ylab="Fst/(1-Fst)", 
     pch=20)

y <- as.numeric(prbi_11_12_comb12_fstlin)
x <- as.numeric(pconnect_avg_comb12_matrix)

mod_pconnect_comb12_matrix <- lm(y~x)
summary(mod_pconnect_comb12_matrix)
#Adjusted R-squared: 0.175, p-value: 0.003409

cor.test(x, y, method='pearson')
#cor: -0.4416733, p-value: 0.003409

plot(pconnect_avg_comb12_matrix, prbi_11_12_comb12_fstlin, 
     xlab="Probability of Larval Dispersal between Populations", ylab="Fst/(1-Fst)", 
     pch=20)
abline(mod_pconnect_comb12_matrix, col="red")

#Partial Mantel test controlling for over water geographic distance

#Load in matrix with over water distances

water_distance_comb12 <- read.csv("OverWater_Comb12_Distance.csv")
water_distance_comb12 <- as.matrix(water_distance_comb12)
water_distance_comb12 <- water_distance_comb12[, -1]
colnames(water_distance_comb12) <- c(1, 7, 8, 9, 10, 11, 19)
rownames(water_distance_comb12) <- c(1, 7, 8, 9, 10, 11, 19)

#Make matrix symmetrical
upperTriangle(water_distance_comb12) <- lowerTriangle(water_distance_comb12, byrow=TRUE)

#Run Partial Mantel test
mantel.partial(pconnect_avg_comb12_matrix, prbi_11_12_comb12_fstlin, water_distance_comb12)
#Mantel statistic r: -0.4151, p-value: 0.963, 5039 permutations



##Pairwise matrix with pops 1 and 2 Fst combined and release site 2480 for pop 1 and pop 19 excluded##

pconnect_avg_comb12_no19_matrix <- pconnect_avg_comb12_matrix[-7, -7]

prbi_11_12_comb12_no19_fstlin <- prbi_11_12_comb12_fstlin[-7, -7]

#Mantel test between average pconnect and linearized Fst 

mantel(pconnect_avg_comb12_no19_matrix, prbi_11_12_comb12_no19_fstlin) 
#Mantel statistic r: -0.3366, p-value: 0.81667, 719 permutations

plot(pconnect_avg_comb12_no19_matrix, prbi_11_12_comb12_no19_fstlin, 
     xlab="Probability of Larval Dispersal between Populations", ylab="Fst/(1-Fst)", 
     pch=20)

y <- as.numeric(prbi_11_12_comb12_no19_fstlin)
x <- as.numeric(pconnect_avg_comb12_no19_matrix)

mod_pconnect_comb12_no19_matrix <- lm(y~x)
summary(mod_pconnect_comb12_no19_matrix)
#Adjusted R-squared: 0.08165, p-value: 0.06892

cor.test(x, y, method='pearson')
#cor: -0.3366305, p-value: 0.06892

plot(pconnect_avg_comb12_no19_matrix, prbi_11_12_comb12_no19_fstlin, 
     xlab="Probability of Larval Dispersal between Populations", ylab="Fst/(1-Fst)", 
     pch=20)
abline(mod_pconnect_comb12_no19_matrix, col="red")

#Partial Mantel test

water_distance_comb12_no19 <- water_distance_comb12[-7, -7]

mantel.partial(pconnect_avg_comb12_no19_matrix, prbi_11_12_comb12_no19_fstlin, water_distance_comb12_no19)

#Mantel statistic r: -0.01532, p-value: 0.52361, 719 permutations




#For reference
colnames(prbi_11_12_fst) <- c(1, 2, 7, 8, 9, 10, 11, 13, 14, 15, 19, 22)
rownames(prbi_11_12_fst) <- c(1, 2, 7, 8, 9, 10, 11, 13, 14, 15, 19, 22)
prbi_geodistance_km[upper.tri(prbi_geodistance_km, diag=T)] = NA

