##Census density estimates for PRBI
#Following PhilippinesAnemonefish2008/Surveys/analyzePhilsurveys 100512.R from PinskyLab GitHub

##Summary stats on surveys

library(tidyverse)

#Read in survey data
surv = read.csv("surveys2009-09-23.density.csv")
glimpse(surv)

# distance
k = surv$RandomSite == 1
k
sum(surv$length[k])/1000 # total length of random site surveys in km (12.44224 km)
sum(surv$length[k])/1000/(252+223)*100 # % of Cebu+Leyte study area 
#(from 1/8/10 calculation in ArcGIS) =2.619%
mean(surv$length[k]) # mean distance of random surveys (=654.85 m)
sd(surv$length[k])/sqrt(sum(k)) #=50.84 m

# area
k = surv$RandomSite == 1
sum(surv$area[k]) # total area in m2 (=111,240.4 m2)
sum(surv$area[k])/1000 #total area in km2 (=111.2404 km2)


###########################################################################

#Density estimates for random sites

##All random sites
k = surv$RandomSite == 1

#Density (fish/km2) of all life stages of PRBI at all random sites using counts/area
(sum(surv$countPRBI[k]))/(sum(surv$area[k]/1000)) #(=0.2517 fish/km2)

#Using dens column of data sheet
mean(surv$densPRBI[k]*1000) #(=0.267 fish/km2)

##Random sites in Cebu
k = surv$Region == "Cebu" & surv$RandomSite == 1

#Density using counts/area
(sum(surv$countPRBI[k]))/(sum(surv$area[k]/1000)) #(=0.20338 fish/km2)

#Density counts using dens column
mean(surv$densPRBI[k]*1000) #(=0.1931 fish/km2)

##Random sites in Leyte
k = surv$Region == "Leyte" & surv$RandomSite == 1

#Density using counts/area
(sum(surv$countPRBI[k]))/(sum(surv$area[k]/1000)) #(=0.3802 fish/km2)

#Density counts using dens column
mean(surv$densPRBI[k]*1000) #(=0.3759 fish/km2)

##Random sites in Danajon
k = surv$Region == "Danajon" & surv$RandomSite == 1

#Density using counts/area
(sum(surv$countPRBI[k]))/(sum(surv$area[k]/1000)) #(=0.2809 fish/km2)

#Density counts using dens column
mean(surv$densPRBI[k]*1000) #(=0.2969 fish/km2)

##############################################################

##Density estimates for adults at random sites

##Adults at all random sites
k = surv$RandomSite == 1

#Density of adults using counts/area
(sum(surv$countPRBIad[k]))/(sum(surv$area[k]/1000)) #(=0.1259 adults/km2)

#Using dens column of data sheet
mean(surv$densPRBIad[k]*1000) #(=0.1289 adults/km2)

##Adults at random sites in Cebu
k = surv$Region == "Cebu" & surv$RandomSite == 1

#Density using counts/area
(sum(surv$countPRBIad[k]))/(sum(surv$area[k]/1000)) #(=0.1017 adults/km2)

#Density counts using dens column
mean(surv$densPRBIad[k]*1000) #(=0.09795 adults/km2)

##Adults at random sites in Leyte
k = surv$Region == "Leyte" & surv$RandomSite == 1

#Density using counts/area
(sum(surv$countPRBIad[k]))/(sum(surv$area[k]/1000)) #(=0.1901 adults/km2)

#Density counts using dens column
mean(surv$densPRBIad[k]*1000) #(=0.18796 adults/km2)

##Adults at random sites in Danajon
k = surv$Region == "Danajon" & surv$RandomSite == 1

#Density using counts/area
(sum(surv$countPRBIad[k]))/(sum(surv$area[k]/1000)) #(=0.1404 adults/km2)

#Density counts using dens column
mean(surv$densPRBIad[k]*1000) #(=0.1137 adults/km2) 

##################################################################

#Checking if each random site has an area value associated with it:
surv_random <- select(surv, RandomSite, area)
surv_random <- filter(surv_random, RandomSite != 0)
surv_random  #Yes, for RandomSite = 1 there is an area at each site

#What about non-random sites?
surv_nonrandom <- select(surv, RandomSite, LinearFishSurvey, MappingFishSurvey, area)
surv_nonrandom <- filter(surv_nonrandom, RandomSite == 0)
#Points with NA for area, have a 0 under RandomSite, Linear FishSurvey, and MappingFishSurvey

#Subset data and filter out the NAs (where there wasn't a survey)
surv_subset <- select(surv, RandomSite, LinearFishSurvey, MappingFishSurvey, Region, length, area, 
                      countPRBI, countPRBIad, densPRBI, densPRBIad)
surv_subset <- filter(surv_subset, area != "NA")


#####################################################################

##Density estimates for all surveys (LinearFishSurvey, MappingFishSurvey, and RandomSites)

##All regions

#Density (fish/km2) of all life stages of PRBI using counts/area
(sum(surv_subset$countPRBI))/(sum(surv_subset$area/1000)) #(=0.4125 fish/km2)

#Using dens column of data sheet
mean(surv_subset$densPRBI*1000) #(=0.9426 fish/km2)

##Cebu
k= surv_subset$Region == "Cebu"

#Density using counts/area
(sum(surv_subset$countPRBI[k]))/(sum(surv_subset$area[k]/1000)) #(=0.3453 fish/km2)

#Density counts using dens column
mean(surv_subset$densPRBI[k]*1000) #(=0.5304 fish/km2)

##Random sites in Leyte
k = surv_subset$Region == "Leyte" 

#Density using counts/area
(sum(surv_subset$countPRBI[k]))/(sum(surv_subset$area[k]/1000)) #(=0.2598 fish/km2)

#Density counts using dens column
mean(surv_subset$densPRBI[k]*1000) #(=0.9152 fish/km2)

##Random sites in Danajon
k = surv_subset$Region == "Danajon" 

#Density using counts/area
(sum(surv_subset$countPRBI[k]))/(sum(surv_subset$area[k]/1000)) #(=1.48 fish/km2)

#Density counts using dens column
mean(surv_subset$densPRBI[k]*1000) #(=3.96 fish/km2)

#################################################################################
##Density estimates for adults in all surveys (LinearFishSurvey, MappingFishSurvey, and RandomSites)

##All regions

#Density (fish/km2) of all life stages of PRBI using counts/area
(sum(surv_subset$countPRBIad))/(sum(surv_subset$area/1000)) #(=0.2510 adults/km2)

#Using dens column of data sheet
mean(surv_subset$densPRBIad*1000) #(=0.6098 adults/km2)

##Cebu
k= surv_subset$Region == "Cebu"

#Density using counts/area
(sum(surv_subset$countPRBIad[k]))/(sum(surv_subset$area[k]/1000)) #(=0.1871 adults/km2)

#Density counts using dens column
mean(surv_subset$densPRBIad[k]*1000) #(=0.2978 adults/km2)

##Random sites in Leyte
k = surv_subset$Region == "Leyte" 

#Density using counts/area
(sum(surv_subset$countPRBIad[k]))/(sum(surv_subset$area[k]/1000)) #(=0.1685 adults/km2)

#Density counts using dens column
mean(surv_subset$densPRBIad[k]*1000) #(=0.6026 adults/km2)

##Random sites in Danajon
k = surv_subset$Region == "Danajon" 

#Density using counts/area
(sum(surv_subset$countPRBIad[k]))/(sum(surv_subset$area[k]/1000)) #(=1.03 adults/km2)

#Density counts using dens column
mean(surv_subset$densPRBIad[k]*1000) #(=2.83 adults/km2)
