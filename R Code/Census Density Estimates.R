##Census density estimates for PRBI
#Following PhilippinesAnemonefish2008/Surveys/analyzePhilsurveys 100512.R from PinskyLab GitHub

##Summary stats on surveys

library(tidyverse)

#Read in survey data from "Downloaded Data" folder
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

IBD_mean *1000000 * 193 / 130 #524 fish/km


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

IBD_lowerCI <- 0.0002 * 1000000 * 193 / 130 #296.9231 fish/km

IBD_upperCI <- 0.0005 * 1000000 * 193 / 130 #742.3077 fish/km


#Density estimate for Cebu (Bohol pops 1 and 2 excluded)

Cebu_dens <- subset(random_surveys, SurveyNum == c(4, 5, 6, 18, 23, 26, 34, 39))

Cebu_dens <- Cebu_dens[-3,]

Cebu_dens$densPRBI #densities at each transect in Cebu (fish/m^2)

hist(Cebu_dens$densPRBI) #doesn't look normal

#Just to compare, ran the mean and se to see what the CI would be if we assumed the distribution is normal
Cebu_mean <- mean(Cebu_dens$densPRBI) #0.0003514099
Cebu_se <- sd(Cebu_dens$densPRBI) / sqrt(length(Cebu_dens$densPRBI)) #9.580003e-05

Cebu_mean *1000000 * 68 / 110 #217 fish/km


##Bootstrapping of Cebu density estimates

library(boot)

x = as.vector(Cebu_dens$densPRBI)

samplemean <- function(x, d) {
  return(mean(x[d]))
}

b = boot(x, samplemean, R=1000)

b
plot(b)

boot.ci(boot.out=b, type="bca") #95% CI 0.0002 - 0.0005

#Convert to square km, then multiply CI by reef area (193km^2) and divide by reef length (130 km)

Cebu_lowerCI <- 0.0002 * 1000000 * 68 / 110 #123.6364 fish/km

Cebu_upperCI <- 0.0005 * 1000000 * 68 / 110 #309.0909 fish/km


