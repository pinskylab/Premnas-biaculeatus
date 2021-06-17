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


#Then plot Fst on y-axis and probability of connection on the x-axis, for now average the two
#potential connectivity values that come up (the matrix isn't symmetrical so there will be 2 values
#between each population)







