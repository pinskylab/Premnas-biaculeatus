## Kernel spreads for other papers
## Generalized gaussian: 
##	PDF = beta/(2*alpha*gamma(1/beta))*exp(-(abs(x)/alpha)^beta)
##		for position x relative to parent
##  standard deviation (= dispersal spread) = sqrt(alpha^3*gamma(3/beta)/gamma(1/beta))

# Catalano et al. 2021 Table S7 Amphiprion clarkii
beta = c(1.03, 0.22, 0.38, 0.67, 5, 0.26, 1.37, 1.49) # years 2012-2018 and then average
alpha = 1/exp(c(-2.36, 4.04, 0.49, -1.52, -3.04, 2.94, -2.32, -2.51))
sqrt(alpha^2*gamma(3/beta)/gamma(1/beta))

# Almany et al. 2017 Amphiprion percula
# parameters from text, Methods Eq. 2 for parameterization
beta = 1
alpha = c(18.9, 13.3) # 2009 and 2011, same as mean dispersal since a Laplacian kernel
sqrt(alpha^2*gamma(3/beta)/gamma(1/beta))

# Almany et al. 2013 Plectropomus areolatus
# Fig. 2 caption for parameters, Supplement p. 11 for parameterization
beta = 3
alpha = (1/4.2e-5)^(1/3)
sqrt(alpha^2*gamma(3/beta)/gamma(1/beta))

# D'Aloia et al. 2015 Elacatinus lori
beta = 1
alpha = 1/0.36
sqrt(alpha^2*gamma(3/beta)/gamma(1/beta))

# Abesamis et al. 2017 Chaetodon vagabundus
# Table 1 for ln(k) and Supplement Eq 3 for parameterization
beta = 3
alpha = (1/exp(-12.8))^(1/3)
sqrt(alpha^2*gamma(3/beta)/gamma(1/beta))

# Almany et al. 2017 Chaetodon vagabundus
# parameters from text, Methods Eq. 2 for parameterization
beta = 1
alpha = 220 # 2011, same as mean dispersal since a Laplacian kernel
sqrt(alpha^2*gamma(3/beta)/gamma(1/beta))

# Williamson et al. 2016 Plectropomus maculatus
beta = 1
alpha = 1/0.00371 # they report 50% retention within 185 km and 95% retention within 811 km. qexp(c(0.5, 0.95), 0.00371) is about right (187, 807)
sqrt(alpha^2*gamma(3/beta)/gamma(1/beta))

# Williamson et al. 2016 Plectropomus leopardus
beta = 1
alpha = 1/0.00628 # they report 50% retention within 110 km and 95% retention within 480 km. qexp(c(0.5, 0.95), 0.00628) is about right (110.4, 477)
sqrt(alpha^2*gamma(3/beta)/gamma(1/beta))

