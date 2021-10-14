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

#Density estimates for random sites (updated 10/7/21)

##All random sites
k = surv$RandomSite == 1

#Sum of survey counts at all random sites
sum(surv$countPRBI[k])  #28 fish

length(surv$countPRBI[k])


##Density estimate at each site where PRBI was found

PRBIdens <- data.frame(density=surv$densPRBI[k])

PRBIdens <- subset(PRBIdens, density > 0)

PRBIdens <- PRBIdens$density * 1000000 #convert from fish/m^2 to fish per km^2

hist(PRBIdens) #doesn't look normal

#Just to compare, ran the mean and se to see what the CI would be if we assumed the distribution is normal
mean(PRBIdens) #634.6966, appears to be skewed by the highest value
sd(PRBIdens) / sqrt(length(PRBIdens)) #std. error is 240.743, makes 95% CI 393.9536-875.4396

mean(PRBIdens) * 300 / 130 #1464.685 fish/km

##Bootstrapping of density estimates

library(boot)

x = as.vector(PRBIdens)

samplemean <- function(x, d) {
  return(mean(x[d]))
}

b = boot(x, samplemean, R=1000)

b
plot(b)

boot.ci(boot.out=b, type="bca") #95% CI 342.6-1352.1

#Multiply CI by reef area (300km^2) and divide by reef length (130 km)

342.6 * 300 / 130 #790.6154 fish/km

1352.1 * 300 / 130 #3120.231 fish/km



##Density estimate for the IBD study region

#Find the surveys in the IBD study region

random_surveys <- subset(surv, surv$RandomSite == 1)
random_surveys

#8 surveys were in our IBD study region: survey numbers 4, 5, 6, 18, 23, 26, 34, 39

IBD_dens <- subset(random_surveys, SurveyNum == c(4, 5, 6, 18, 23, 26, 34, 39))

IBD_dens$densPRBI #densities at each transect in IBD study region (fish/m^2)

hist(IBD_dens$densPRBI) #doesn't look normal

#Just to compare, ran the mean and se to see what the CI would be if we assumed the distribution is normal
IBD_mean <- mean(IBD_dens$densPRBI) #0.0003527602
IBD_se <- sd(IBD_dens$densPRBI) / sqrt(length(IBD_dens$densPRBI)) #8.297625e-05

IBD_mean *1000000 * 300 / 130 #814 fish/km


##Bootstrapping of IBD density estimates

library(boot)

x = as.vector(IBD_dens$densPRBI)

samplemean <- function(x, d) {
  return(mean(x[d]))
}

b = boot(x, samplemean, R=1000)

b
plot(b)

boot.ci(boot.out=b, type="bca") #95% CI 0.0002 - 0.0005

#Convert to square km, then multiply CI by reef area (300km^2) and divide by reef length (130 km)

IBD_lowerCI <- 0.0002 * 1000000 * 300 / 130 #461.5385 fish/km

IBD_upperCI <- 0.0005 * 1000000 * 300 / 130 #1153.846 fish/km



#####################################################################



##Old code

#Confidence interval for number of fish
hist(surv$countPRBI[k]) #does not look normal, will bootstrap instead

random_density <- tibble(count=surv$countPRBI[k])
random_density

install.packages("boot")
library(boot)

library(tidyverse) 
library(infer)

surv_resample <- rep_sample_n(random_density, size = 1, replace = TRUE, reps = 10000)

surv_resample
mean(surv_resample$count) #1.4546

#Compare bootstrap mean to sample mean
mean(surv$countPRBI[k]) #1.473684

#Generate confidence interval:
#Calculate standard error
std_error <- sd(surv_resample$count) / sqrt(length(surv_resample$count))
std_error #0.02213285

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
