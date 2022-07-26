######################################################################
## Resampling approach to estimate 95% CI on dispersal kernel spread
######################################################################

# Read in scripts (found in "Scripts for Sigma Estimates" folder)
source('DeFromNb.R')
source('sigmaFrom_m.r')
source('summarizeSigmas.R')

# Set parameters
niter <- 1000000 # number of iterations
alpha <- 1.678 # minimum age of maturity (for Waples Nb to Ne calculations)
AL <- 7.052 # adult lifespan (for Waples Nb to Ne calculations)
A <- 130 # length of reef (km)
m <- 6.253e-05 # isolation by distance slope (genetic distance/km)
mse <- 2.518e-05 # standard error on m


# Estimate sigma (dispersal kernel spread)

#With over-water distance, df=3 and Waples_Nbu95=69,421:

# Set parameters
niter <- 1000000 # number of iterations
alpha <- 1.678 # minimum age of maturity (for Waples Nb to Ne calculations)
AL <- 7.052 # adult lifespan (for Waples Nb to Ne calculations)
A <- 130 # length of reef (km)
m <- 5.393e-05 # isolation by distance slope (genetic distance/km)
mse <- 1.663e-05 # standard error on m

D1 <- DeFromNb(niter, AL, alpha, Waples_Nb=6942.1, Waples_Nbl95=779.5, Waples_Nbu95=69421, A=130) # estimate De (58.6325)
s1 <- sigmaFrom_m(De=D1$De, Des=D1$Des, m=m, mse=mse, dims=1) # estimate sigma
s1$sigma_point # report point value for sigma (8.891716)
summarizeSigmas(sigmas=s1$sigmas, sg=3) # report distribution for sigma

#Point value for sigma: 8.891716, mean: 8.55, median: 7.92, SD: 4.79, Lower 95% CI: 2.3, Upper 95% CI: 18.4

#Calculate number of bootstrap values above 11km (A.clarkii dispersal spread estimate)
s1_df <- as.data.frame(s1)
s1_11 <- subset(s1_df, s1_df$sigmas>11)
s1_11
length(s1_11$sigmas) #234,583 out of 1,000,000 iterations are greater than 11 km
summary(s1_11)
234583/1000000 #0.2346

#Calculate number of bootstrap values above 17km (A. percula dispersal spread estimate)
s1_17 <- subset(s1_df, s1_df$sigmas>17)
s1_17
length(s1_17$sigmas) #37,857 out of 1,000,000 iterations are greater than 17 km
summary(s1_17)
37857/1000000 #0.037857

###########################################################

#Estimating adult lifespan (AL) and age at maturity (alpha)

#Using FishLife package

library(devtools)
devtools::install_github("james-thorson/FishLife")
library(FishLife) #version 2.0.0

vignette("tutorial", "FishLife")

#Description of life history parameters in FishLife model:
#Lmax- maximum length
#Linf (looks like Loo on graph)- asymptotic maximum length
#Winfinity- asymptotic maximum mass
#K- Brody growth coefficient
#Lmat- length at maturity (Lm in R package?)
#Amat- age at 50% sexual maturity, in years (tm in R package?)
#Amax- maximum age, in years (tmax in R package?)
#M- instantaneous natural mortality rate
#Temp- average environmental temperature (degrees Celsius)
#t0- age at which length would be zero according to estimated von Bertalanffy growth curve

#Adult lifespan (AL) = Amax
#Age at maturity (alpha) = Amat

#Generate basic plot for Premnas biaculeatus

Predict = Plot_taxa(Search_species(Genus="Premnas",Species="biaculeatus")$match_taxonomy, mfrow=c(2,2))
#Not generating plots correctly

#Predict means for life history parameters (in log-space except for temperature)
knitr::kable(Predict[[1]]$Mean_pred, digits=3)

#Convert to predictive median by exponentiating the variables:
knitr::kable(c(exp(Predict[[1]]$Mean_pred[-8]),Predict[[1]]$Mean_pred['Temperature']), digits=3)
#tmax is 5.598, tm is 1.297

#Calculate predictive mean by exponentiating and bias-correcting:
knitr::kable(c(exp(Predict[[1]]$Mean_pred[-8]+0.5*diag(Predict[[1]]$Cov_pred)[-8]),Predict[[1]]$Mean_pred['Temperature']), digits=3)
#tmax is 7.052, tm is 1.678

################################################################
#Correcting Ne estimate for overlapping generations
#Using equations from Waples et al. 2014

#1st use Nb/Ne = 0.485 + 0.758*log(AL/alpha)

Nb_Ne <- 0.485 + (0.758 * log10(7.052/1.678))
Nb_Ne  #0.958

#2nd calculate adjusted Nb using Nbadj = Nb / (1.26 - 0.323*(Nb/Ne))
Nb <- 6942.1
Nb_adj <- Nb / (1.26 - 0.323 * Nb_Ne)
Nb_adj  #7302.201

#3rd calculate Ne using Ne=Nbadj / Nb/Ne

Ne <- Nb_adj / Nb_Ne
Ne #7625.296
Ne/130 #58.65613


##Correcting lower and upper bounds of confidence interval

Nbl <- 779.5
Nb_adj <- Nbl / (1.26 - 0.323 * Nb_Ne)
Nb_adj #819.9342

Nel <- Nb_adj / Nb_Ne
Nel #856.2133
Nel/130 #6.586256

Nbu <- 69421
Nb_adj <- Nbu / (1.26 - 0.323 * Nb_Ne)
Nb_adj #73022.01

Neu <- Nb_adj / Nb_Ne
Neu #76252.96
Neu/130 #586.5613
