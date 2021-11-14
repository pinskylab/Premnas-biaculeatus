######################################################################
## Resampling approach to estimate 95% CI on dispersal kernel spread
######################################################################

# Read in scripts
source('DeFromD.R')
source('DeFromNb.R')
source('sigmaFrom_m.r')
source('sigmaFrom_NbWright.r')
source('summarizeSigmas.R')

# Set parameters
niter <- 1000000 # number of iterations
alpha <- 1.678 # minimum age of maturity (for Waples Nb to Ne calculations)
AL <- 7.052 # adult lifespan (for Waples Nb to Ne calculations)
A <- 130 # length of reef (km)
m <- 6.253e-05 # isolation by distance slope (genetic distance/km)
mse <- 2.518e-05 # standard error on m


# Estimate sigma (dispersal kernel spread) in various ways

	# Using Waples Nb (# breeders) and IBD slope
	# Waples_Nb value and JackKnife CI from NeEstimator results (based on lowest allele frequency of 0.02)
D1 <- DeFromNb(niter, AL, alpha, Waples_Nb=6942.1, Waples_Nbl95=779.5, Waples_Nbu95=10000, A=130) # estimate De
s1 <- sigmaFrom_m(De=D1$De, Des=D1$Des, m=m, mse=mse, dims=1) # estimate sigma
s1$sigma_point # report point value for sigma (8.2577) with df=129
summarizeSigmas(sigmas=s1$sigmas, sg=3) # report distribution for sigma
#With df=129:
#Mean- 9.01, median- 8.23, SD- 6.11, Lower 95% CI- 5.48, Upper 95% CI- 16.8

#With df=3:
#Point value for sigma: 8.257652, mean: 9.02, median: 8.23, SD: 8.32, Lower 95% CI: 5.48, Upper 95% CI: 16.8

#With df=3 and Waples_Nbu95=69,421:

D1 <- DeFromNb(niter, AL, alpha, Waples_Nb=6942.1, Waples_Nbl95=779.5, Waples_Nbu95=69421, A=130) # estimate De
s1 <- sigmaFrom_m(De=D1$De, Des=D1$Des, m=m, mse=mse, dims=1) # estimate sigma
s1$sigma_point # report point value for sigma (8.2577)
summarizeSigmas(sigmas=s1$sigmas, sg=3) # report distribution for sigma

#Point value for sigma: 8.257652, mean: 8.33, median: 7.41, SD: 6.36, Lower 95% CI: 2.11, Upper 95% CI: 19.4


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

#95% CI for De
Del95 = quantile(Des, 0.025) #49.75894
Deu95 = quantile(Des, 0.975) #74.81362

#Mean and se for De

de_mean <- mean(Des) #61.59
de_se <- sd(Des) / sqrt(length(Des)) #2.835

de_mean + 1.96*de_se #67.15
de_mean - 1.96*de_se #56.04


	# Using Waples Nb and Migraine estimate of Wright's neighborhood size
s2 <- sigmaFrom_NbWright(De=D1$De, Des=D1$Des, NbWright=7738, NbWright_l95=6806, NbWright_u95=10527, dims=1)
s2$sigma_point
summarizeSigmas(sigmas=s2$sigmas, sg=3)
			
	# Erroneous calculation using census density   #Not sure what se would be on census density estimate
D3 <- DeFromD(niter, D=188.1, Dse=45, nen=1, nense=0)
s3 <- sigmaFrom_m(De=D3$De, Des=D3$Des, m=m, mse=mse, dims=1)
s3$sigma_point
summarizeSigmas(sigmas=s3$sigmas, sg=3)
	
	# Erroneous calculation using 0.1% of census density
D4 <- DeFromD(niter, D=188.1, Dse=45, nen=0.001, nense=0)
s4 <- sigmaFrom_m(De=D4$De, Des=D4$Des, m=m, mse=mse, dims=1)
s4$sigma_point
summarizeSigmas(sigmas=s4$sigmas, sg=3)


###########################################################

#Estimating adult lifespan (AL) and age at maturity (alpha)

#Using FishLife package

library(devtools)
devtools::install_github("james-thorson/FishLife")
library(FishLife)

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
