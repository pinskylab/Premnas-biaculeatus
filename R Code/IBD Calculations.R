##IBD Calculations


##2009-11-12 data without populations 13, 14, 15, 19, 22 (East-West transect wout 19)

#dispersal= sqrt(1/4Dem)

#Slope= 6.253e-05 from linear regression
#So m=6.253e-05 for the IBD equation

m <- 6.253e-05

## Using Ne from Ne estimator: lowest allele frequency used is 0.02
#Estimated Ne=6942.1, approximate as 6942

#Area of coral reef habitat: 300km2

#De=Ne/Area
#De=6942/300
de <- 6942/300  #De=23.14 fish/km

disp <- sqrt(1/(4*de*m))  #13.1445 km


#Multiplying by the length of the reef 
#Length of reef (covering locations of Pops 1-2, 7-11): 130 km

de_lin <- de*130

disp_lin <- sqrt(1/(4*de_lin*m))  #1.15 km


#Dividing Ne by the length of the reef

de_length <- 6942/130  #De=53.4 fish/km

disp_length <- sqrt(1/(4*de_length*m))  #8.653 km 


##############################################################

##Using census density estimate (all PRBI at random sites)

de_cen <- 0.2517  #Census density=0.2517 fish/km

disp_cen <- sqrt(1/(4*de_cen*m))  #126.033 km


#Multiplying by the length of the reef 

de_cen_lin <- de_cen*130

disp_cen_lin <- sqrt(1/(4*de_cen_lin*m))  #11.054 km 


#Dividing census sum of fish counts by length of reef

#28 fish counted at random sites

de_cen_length <- 28/130  #0.215 fish/km

disp_cen_length <- sqrt(1/(4*de_cen_length*m)) #136.24km

