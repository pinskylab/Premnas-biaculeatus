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
(sum(surv$countPRBIad[k]))/(sum(surv$area[k]/1000)) #(=0.1404 adults/km2)

#Using dens column of data sheet
mean(surv$densPRBIad[k]*1000) #(=0.1137 adults/km2)

##Adults at random sites in Cebu
k = surv$Region == "Cebu" & surv$RandomSite == 1

#Density using counts/area
(sum(surv$countPRBIad[k]))/(sum(surv$area[k]/1000)) #(=0.1017 fish/km2)

#Density counts using dens column
mean(surv$densPRBIad[k]*1000) #(=0.09795 fish/km2)

##Adults at random sites in Leyte
k = surv$Region == "Leyte" & surv$RandomSite == 1

#Density using counts/area
(sum(surv$countPRBIad[k]))/(sum(surv$area[k]/1000)) #(=0.1901 fish/km2)

#Density counts using dens column
mean(surv$densPRBIad[k]*1000) #(=0.18796 fish/km2)

##Adults at random sites in Danajon
k = surv$Region == "Danajon" & surv$RandomSite == 1

#Density using counts/area
(sum(surv$countPRBIad[k]))/(sum(surv$area[k]/1000)) #(=0.1404 fish/km2)

#Density counts using dens column
mean(surv$densPRBIad[k]*1000) #(=0.1137 fish/km2) 


##Random sites: Fish/km for PRBI, by latitude, assuming 150m wide reef

par(mfrow=c(1,3))
par(cex=.7, cex.lab = 1.4)
ylims = c(0, 400)

#Cebu
k = surv$Region == "Cebu" & surv$RandomSite == 1
plot(surv$lat[k], surv$densPRBI[k]*150*1000, main = "P. biaculeatus in Cebu\n(random sites, 150m reef)", 
     xlab = "Latitude (°)", ylab = "fish per km", ylim= ylims)
abline(h = mean(surv$densPRBI[k]*150*1000), lty=3)
printMeanSE(surv$densPRBI[k]*150*1000)

#Leyte
k = surv$Region == "Leyte" & surv$RandomSite == 1
plot(surv$lat[k], surv$densPRBI[k]*150*1000, main = "P. biaculeatus in Leyte\n(random sites, 150m reef)", 
     xlab = "Latitude (°)", ylab = "fish per km", ylim= ylims)
abline(h = mean(surv$densPRBI[k]*150*1000), lty=3)
printMeanSE(surv$densPRBI[k]*150*1000)

#Danajon
k = surv$Region == "Danajon" & surv$RandomSite == 1
plot(surv$lat[k], surv$densPRBI[k]*150*1000, 
     main = "P. biaculeatus on Danajon Bank\n(random sites, 150m reef)", 
     xlab = "Latitude (°)", ylab = "fish per km", ylim= ylims)
abline(h = mean(surv$densPRBI[k]*150*1000), lty=3)
printMeanSE(surv$densPRBI[k]*150*1000)