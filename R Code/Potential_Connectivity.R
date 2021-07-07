## Potentially Connectivity Matrices
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

latlon_matrix <- nc_open("coral_all_release_site_coarse.nc")
print(latlon_matrix)

lon_data <- ncvar_get(latlon_matrix, varid="lon")
print(lon_data)
dim(lon_data)
lon_data[1:10, 1:10]

lat_data <- ncvar_get(latlon_matrix, varid="lat")
print(lat_data)

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

i <- which(c(release_site_data)==2427)
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
#Pop 11- site 2555 (maybe)
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
mod_pconnect_avg <- lm(pop19$GeneticDistance ~ pconnect_avg_df$pconnect_avg)
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
summary(mod_pconnect_dst) #Adjusted R-squared: -0.1227, p:0.583

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
  geom_smooth(method="lm", se=FALSE, color="steelblue") +
  theme_bw() +
  xlab("PConnect- Pop 19 Destination") +
  ylab("Fst/(1-Fst)")

#Could plot with a legend showing the population number 
#Could map coordinates used from the pconnect matrix and compare to our populations, 
#add arrows to map weighted based on pconnect




