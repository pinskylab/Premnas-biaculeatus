##PRBI IBD Plots##

library(tidyverse)
library(broom)
library(wesanderson)
library(formattable)

#IBD Pattern Plot

#Convert genetic distance and geographic distance matrices into a dataframe

distanceframe <- tibble(GeographicDistance=as.vector(prbi_EW_no19_geodist_km), 
                        GeneticDistance=as.vector(prbi_11_12_EW_no19_fstlin))
distanceframe <- drop_na(distanceframe)

mod_distanceframe <- lm(distanceframe$GeneticDistance~distanceframe$GeographicDistance)
summary(mod_distanceframe)

#Calculating 95% CI on slope
coef(mod_distanceframe)
slope <- coef(mod_distanceframe)[2] #slope

summ <- summary(mod_distanceframe)
se <- summ$coefficients[2,2] #standard error on slope estimate

upperci <- slope + 1.96*se
lowerci <- slope - 1.96*se

confint(mod_distanceframe)

#Using confidence intervals from augment function
augmented_mod_distanceframe <- augment(mod_distanceframe, interval = 'confidence')
glimpse(augmented_mod_distanceframe)

ggplot(data=augmented_mod_distanceframe, aes(x=distanceframe$GeographicDistance, 
                                             y=distanceframe$GeneticDistance)) +
  geom_point(color="slateblue", size=1.75) +
  geom_line(aes(y=.fitted), color="slateblue", size=1.5) +
  geom_ribbon(aes(ymin=.lower, ymax=.upper), alpha=0.3, fill="slateblue") +
  theme_bw() +
  xlab("Geographic Distance (km)") +
  ylab("Genetic Distance (Fst/(1-Fst)")


ggplot(data=distanceframe, aes(x=GeographicDistance, y=GeneticDistance)) +
  geom_point()


#Data Tables

Ne_sigmaestimate <- data.frame(Metric=c("Slope", "Ne", "Reef Length", "De", "Sigma Estimate"),
                          Value=c(6.253e-05, 6942, 130, 53.4, 8.653))
Ne_sigmaestimate
