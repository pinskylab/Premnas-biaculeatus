## Potentially Connectivity Matrices
#From Diane Thompson- for 10-day PLD

##Loading in matrices with tidync package

install.packages("tidync")
library(tidync)

library(ncdf4)
connect_matrix <- nc_open("CORAL-connect_25_historical_pconnect_10-day_noedges2.nc")
print(connect_matrix)  #This gives the meta-data for the matrix
#[dst_site, src_site] - dst=destination?, src=?

connect_data <- ncvar_get(connect_matrix, varid="pconnect") #this is the actual potential connectivity matrix
print(connect_data)

latlon_matrix <- nc_open("coral_all_release_site_coarse.nc")
print(latlon_matrix)

lon_data <- ncvar_get(latlon_matrix, varid="lon")
print(lon_data)

lat_data <- ncvar_get(latlon_matrix, varid="lat")
print(lat_data)

release_site_data <- ncvar_get(latlon_matrix, varid="release_site")
print(release_site_data)

#Set path and filename
ncpath <- "C:\Users\kyras\OneDrive\Documents\GitHub\Premnas-biaculeatus\Downloaded Data\"





