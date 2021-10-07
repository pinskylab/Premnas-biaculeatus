##PRBI IBD Plots##

library(tidyverse)
library(broom)
library(wesanderson)
library(formattable)

#IBD Pattern Plot

#Convert genetic distance and geographic distance matrices into a dataframe

distanceframe <- tibble(GeographicDistance=as.vector(prbi_EW_no19_geodist_km), 
                        GeneticDistance=as.vector(prbi_11_12_EW_no19_fstlin), 
                        WaterDistance=as.vector(water_distance_no19))
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

#Plot showing populations vs. pop 19

pop19 <- tibble(GeographicDistance = as.vector(water_distance[8,]),
                GeneticDistance = as.vector(prbi_11_12_EW_fstlin[8,]))
pop19 <- drop_na(pop19)

ggplot(data=distanceframe, aes(x=distanceframe$GeographicDistance, 
                                             y=distanceframe$GeneticDistance)) +
  geom_point(color="#7b3294", size=1.75) +
  geom_point(data=pop19, aes(x=GeographicDistance, y=GeneticDistance), color="#008837", size=1.75) +
  geom_smooth(method="lm", color="#7b3294", fill="#c2a5cf") +
  theme_bw() +
  xlab("Geographic Distance (km)") +
  ylab("Genetic Distance (Fst/(1-Fst)")

#Data Tables

Ne_sigmaestimate <- data.frame(Metric=c("Slope", "Ne", "Reef Length", "De", "Sigma Estimate"),
                          Value=c(6.253e-05, 6942, 130, 53.4, 8.653))
Ne_sigmaestimate



##Plot for comparing dispersal distances across coral reef fish species##

#Load in data

fish_dispersal <- read_csv("Fish_DispersalDistances.csv")
fish_dispersal <- rename(fish_dispersal, Dispersal_Distance = "Dispersal Distance (km)")
fish_dispersal$Species[7] <- "P. maculatus"

ggplot(fish_dispersal, aes(x=Species, y=Dispersal_Distance)) +
  geom_point( color=Species, size=4, alpha=0.6)

fish_dispersal %>%
  mutate(Species = fct_relevel(Species, 
                               "P. leopardus", "P. maculatus", "C. vagabundus", "A. percula", "P. areolatus", 
                               "A. clarkii", "P. biaculeatus","E. lori")) %>%
ggplot(aes(x=Species, y=Dispersal_Distance)) +
  geom_segment( aes(x=Species, xend=Species, y=0, yend=Dispersal_Distance, color=Species), size=2) +
  geom_point(aes(color=Species), size=4, alpha=0.6) +
  theme_light() +
  coord_flip() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position="none"
  ) +
  ylab("Dispersal Distance (km)")

#Graph without P. leopardus and P. areolatus

fish_dispersal_short <- fish_dispersal[-c(4, 8), ]

fish_dispersal_short %>%
  mutate(Species = fct_relevel(Species, 
                              "P. maculatus", "C. vagabundus", "A. percula",
                               "A. clarkii", "P. biaculeatus","E. lori")) %>%
  ggplot(aes(x=Species, y=Dispersal_Distance)) +
  geom_segment( aes(x=Species, xend=Species, y=0, yend=Dispersal_Distance, color=Species), size=2) +
  geom_point(aes(color=Species), size=4, alpha=0.6) +
  theme_light() +
  coord_flip() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position="none"
  ) +
  ylab("Dispersal Distance (km)")
