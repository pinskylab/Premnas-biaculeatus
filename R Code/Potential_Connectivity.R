## Potential Connectivity Matrices
#From Diane Thompson- for 10-day PLD

##Loading in matrices from "Downloaded_Data" folder

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

#Indexing to identify release sites that match the coordinates of our sampling sites
c(lat_data)[i]  #lat of release site 1
c(lon_data)[i]  #lon of release site 1

plot(c(lon_data)[i], c(lat_data)[i]) 

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


##Pairwise matrix with pops 1 and 2 Fst combined and release site 2480 for pop 1#
#Combining pops 1 and 2 as they are about 10 km apart, meaning they are too close together
#for two distinct release sites

#Pop 1- site 2480
#Pop 7- site 2401
#Pop 8- site 2426
#Pop 9- site 2427
#Pop 10- site 2454 
#Pop 11- site 2455 
#Pop 19- site 2483

#Fst matrix is read in and linearized in "PRBI IBD analysis.R" script

pconnect_avg_comb12_matrix <- pconnect_avg_matrix[-2, -2]

upperTriangle(prbi_11_12_comb12_fstlin) <- lowerTriangle(prbi_11_12_comb12_fstlin, byrow=TRUE)

#Mantel test between average pconnect and linearized Fst 

library(vegan)
mantel(pconnect_avg_comb12_matrix, prbi_11_12_comb12_fstlin) 
#Mantel statistic r: -0.4417, p-value: 0.975, 5039 permutations

#Get 95% CI
cor.test(pconnect_avg_comb12_matrix[lower.tri(pconnect_avg_comb12_matrix)], 
         prbi_11_12_comb12_fstlin[lower.tri(prbi_11_12_comb12_fstlin)])
#-0.734 - -0.012

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

#Load in matrix with over water distances from "Coordinates_Distances" folder

water_distance_comb12 <- read.csv("OverWater_Comb12_Distance.csv")
water_distance_comb12 <- as.matrix(water_distance_comb12)
water_distance_comb12 <- water_distance_comb12[, -1]
colnames(water_distance_comb12) <- c(1, 7, 8, 9, 10, 11, 19)
rownames(water_distance_comb12) <- c(1, 7, 8, 9, 10, 11, 19)

#Make matrix symmetrical
upperTriangle(water_distance_comb12) <- lowerTriangle(water_distance_comb12, byrow=TRUE)

#Run Partial Mantel test with geographic distance as control
mantel.partial(pconnect_avg_comb12_matrix, prbi_11_12_comb12_fstlin, water_distance_comb12)
#Mantel statistic r: -0.4151, p-value: 0.963, 5039 permutations

#Run Partial Mantel test with pconnect as control
mantel.partial(water_distance_comb12, prbi_11_12_comb12_fstlin, pconnect_avg_comb12_matrix)
#Mantel statistic r: -0.06837, p-value: 0.591, 5039 permutations


#Mantel test excluding pops 8, 9, 10

pconnect_avg_comb12_no8_9_10_matrix <- pconnect_avg_comb12_matrix[c(-3, -4, -5), c(-3, -4, -5)]

prbi_11_12_comb12_no8_9_10_fstlin <- prbi_11_12_comb12_fstlin[c(-3, -4, -5), c(-3, -4, -5)]

mantel(pconnect_avg_comb12_no8_9_10_matrix, prbi_11_12_comb12_no8_9_10_fstlin, permutations=999)
#R: -0.6616, p: 0.91667, 23 permutations

#Get 95% CI
cor.test(pconnect_avg_comb12_no8_9_10_matrix[lower.tri(pconnect_avg_comb12_no8_9_10_matrix)], 
         prbi_11_12_comb12_no8_9_10_fstlin[lower.tri(prbi_11_12_comb12_no8_9_10_fstlin)])
#-0.959 - 0.324

plot(pconnect_avg_comb12_no8_9_10_matrix, prbi_11_12_comb12_no8_9_10_fstlin, 
     xlab="Probability of Larval Dispersal between Populations", ylab="Fst/(1-Fst)", 
     pch=20)


##Pairwise matrix with pops 1 and 2 Fst combined and release site 2480 for pop 1 and pop 19 excluded##

pconnect_avg_comb12_no19_matrix <- pconnect_avg_comb12_matrix[-7, -7]

prbi_11_12_comb12_no19_fstlin <- prbi_11_12_comb12_fstlin[-7, -7]

#Mantel test between average pconnect and linearized Fst 

mantel(pconnect_avg_comb12_no19_matrix, prbi_11_12_comb12_no19_fstlin) 
#Mantel statistic r: -0.3366, p-value: 0.81667, 719 permutations

cor.test(pconnect_avg_comb12_no19_matrix[lower.tri(pconnect_avg_comb12_no19_matrix)], 
         prbi_11_12_comb12_no19_fstlin[lower.tri(prbi_11_12_comb12_no19_fstlin)])
#-0.724 - 0.212

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

#Partial Mantel test controlling for geographic distance

water_distance_comb12_no19 <- water_distance_comb12[-7, -7]

mantel.partial(pconnect_avg_comb12_no19_matrix, prbi_11_12_comb12_no19_fstlin, water_distance_comb12_no19)

#Mantel statistic r: -0.01532, p-value: 0.52361, 719 permutations

#Partial Mantel controlling for pconnect

mantel.partial(water_distance_comb12_no19, prbi_11_12_comb12_no19_fstlin, pconnect_avg_comb12_no19_matrix)
#Mantel statistic r: 0.4026, p-value: 0.098611, 719 permutations

#Plot showing 19 in different color

pconnect_avg_comb12_matrix[upper.tri(pconnect_avg_comb12_matrix, diag=T)] = NA
prbi_11_12_comb12_fstlin[upper.tri(prbi_11_12_comb12_fstlin, diag=T)] = NA


pconnectframe_comb12 <- tibble(PConnect=as.vector(pconnect_avg_comb12_no19_matrix), 
                        GeneticDistance=as.vector(prbi_11_12_comb12_no19_fstlin))
pconnectframe_comb12 <- drop_na(pconnectframe_comb12)

pop19_comb12 <- tibble(PConnect = as.vector(pconnect_avg_comb12_matrix[7,]),
                GeneticDistance = as.vector(prbi_11_12_comb12_fstlin[7,]))
pop19_comb12 <- drop_na(pop19_comb12)

#Separate symbol for 8&9, 8&10, 9&10

pop8910 <- tibble(PConnect= as.vector(c(pconnect_avg_comb12_matrix[4, 3], pconnect_avg_comb12_matrix[5, 3],
                                        pconnect_avg_comb12_matrix[5, 4])),
                  GeneticDistance= as.vector(c(prbi_11_12_comb12_fstlin[4, 3], prbi_11_12_comb12_fstlin[5, 3],
                                               prbi_11_12_comb12_fstlin[5, 4])))


ggplot(data=pconnectframe_comb12, aes(x=PConnect, 
                                             y=GeneticDistance)) +
  geom_point(size=1.75, color="#7b3294") +
  geom_point(data=pop19_comb12, aes(x=PConnect, y=GeneticDistance), color="#008837", size=1.75) +
  geom_smooth(data=pop19_comb12, method="lm", se=FALSE, color="#008837") +
  theme_bw() +
  xlab("Potential Connectivity") +
  ylab("Genetic Distance (Fst/(1-Fst)")

#With points and circles

ggplot(data=pconnectframe_comb12, aes(x=PConnect, 
                                      y=GeneticDistance)) +
  geom_point(size=2.5) +
  geom_point(data=pop19_comb12, aes(x=PConnect, y=GeneticDistance), size=3, shape=1) +
  geom_smooth(data=pop19_comb12, method="lm", se=FALSE, color="black", linetype="dotted", size=1) +
  geom_point(data=pop8910, aes(x=PConnect, y=GeneticDistance), size=3, shape=15) +
  theme_bw() +
  xlab("Potential Connectivity") +
  ylab("Genetic Distance (Fst/(1-Fst)") +
  annotate("text", x = 0.038, y = -0.005, label = "8") +
  annotate("text", x = 0.023, y = -0.0023, label = "7") +
  annotate("text", x = 0.0062, y = 0.0023, label = "9") +
  annotate("text", x = 0.0045, y = 0.011, label = "1&2") +
  annotate("text", x = -0.0045, y = 0.0075, label = "10") +
  annotate("text", x = 0.007, y = 0.009, label = "11") +
  theme(axis.text = element_text(size = 10))


##Plot Pop 19 compared to other pops with pops 1 and 2 combined (release site 2480)##

#Make vector with pconnect between pop 19 and other pops

#Pop 1- site 2480
#Pop 7- site 2401
#Pop 8- site 2426
#Pop 9- site 2427
#Pop 10- site 2454 
#Pop 11- site 2455 
#Pop 19- site 2483

#1st in direction of pop 19 to other pops
#[dst site, src site] would be [other pop, pop 19], so [other site, 2483]
#Order of vector: Pop 1, 7, 8, 9, 10, 11, 19 
pconnect_pop19src_comb12 <- c(connect_data[2480, 2483],connect_data[2401, 2483], connect_data[2426, 2483],
                            connect_data[2427, 2483], connect_data[2454, 2483], connect_data[2455, 2483])
pconnect_pop19src_comb12

#2nd in direction of other pops to pop 19
pconnect_pop19dst_comb12 <- c(connect_data[2483, 2480], connect_data[2483, 2401], connect_data[2483, 2426],
                            connect_data[2483, 2427], connect_data[2483, 2454], connect_data[2483, 2455])

pconnect_pop19dst_comb12

pconnect_avg_comb12 <- apply(rbind(pconnect_pop19src_comb12, pconnect_pop19dst_comb12), 2, mean)
pconnect_avg_comb12

pop_numbers_comb12 <- c(1, 7, 8, 9, 10, 11)

pop19_comb12_linfst <- as.vector(prbi_11_12_comb12_fstlin[7,])
pop19_comb12_linfst <- pop19_comb12_linfst[-7]

        
pconnect_avg_comb12_df <- data.frame(pop_numbers_comb12, pconnect_avg_comb12, pop19_comb12_linfst)
pconnect_avg_comb12_df


#Plotting average pconnect vs. lin Fst
mod_pconnect_comb12_avg <- lm(pconnect_avg_comb12_df$pop19_comb12_linfst ~ pconnect_avg_comb12_df$pconnect_avg_comb12)
summary(mod_pconnect_comb12_avg) #Adjusted R-squared:0.900235922, p:0.0024

mantel(pconnect_avg_comb12_df$pconnect_avg_comb12, pconnect_avg_comb12_df$pop19_comb12_linfst)

cor.test(pconnect_avg_comb12_df$pconnect_avg_comb12, pconnect_avg_comb12_df$pop19_comb12_linfst, method = 'pearson')
#cor: -0.9600763, p-value: 0.002359, df=4

plot(pconnect_avg_comb12_df$pconnect_avg_comb12, pconnect_avg_comb12_df$pop19_comb12_linfst, pch=19)
abline(mod_pconnect_comb12_avg, col="red")

pconnect_avg_comb12_df$pop_numbers_comb12 <- as.factor(pconnect_avg_comb12_df$pop_numbers_comb12)
class(pconnect_avg_comb12_df$pop_numbers_comb12)

ggplot(data=pconnect_avg_comb12_df, aes(x=pconnect_avg_comb12, 
                                      y=pop19_comb12_linfst, color=pop_numbers_comb12)) +
  geom_point(size=2.5) +
  geom_smooth(method="lm", se=FALSE, color="steelblue") +
  theme_bw() +
  xlab("PConnect") +
  ylab("Fst/(1-Fst)")


#Dispersal spread estimate using equation from Siegel et al. 2003

#standard deviation of current velocity in IBD region: 0.2004 m/s, need to convert to km/day

sd_currentvel <- (0.2004 * 86400) / 1000 #17.31456 km/day


pconnect_spread <- 2.238 * sd_currentvel * 10^0.5 #122.5382

(2.238 * sd_currentvel * 10)^0.5

#But if the power to 0.5 is applied after, answer is 19.685 km


