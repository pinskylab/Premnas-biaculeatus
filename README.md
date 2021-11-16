# Premnas biaculeatus
 
Data files and analysis for the manuscript "Isolation by distance and isolation by oceanography in Maroon Anemonefish (Premnas biaculeatus)" focused on estimating the larval dispersal kernels of Premnas biaculeatus using Isolation-by-Distance genetic estimates and comparing genetic differentiation to probability of dispersal based on oceanographic simulations.

The repository is organized in the following folders:

Downloaded_Data: raw data files, unmanipulated by the authors.
   "Amphiprion_Microsats.csv" details on microsatellites screened for use. Used to develop Tables 1 and S1.
    "coral_all_release_site_coarse.nc" from Thompson et al. 2018. Latitude and longitude for particle release sites used in potential connectivity analyses.
    "CORAL-connect_25_historical_pconnect_10-day_noedges2.nc" from Thompson et al. 2018. Potential connectivity between particle release sites.
     "PRBI_GenAlEx_2009-11-12.csv" genotypes for all Premnas biaculeatus individuals at 16 microsatellite loci in GenAlEx format.
     "PRBI_genepop_2009_11-12.gen.txt" genotypes for all Premnas biaculeatus individuals at 16 microsatellite loci in Genepop format. Used as the input file for Fst                    calculations.
     "surveys2009-09-23.density.csv" number of anemonefish, including Premnas biaculeatus, found at sampling sites in the central Philippines. Used for census density                 calculations.
    
 
Data: data files that have been formatted or manipulated by authors, in addition to output files from running Genepop functions in RStudio. Contains 3 subfolders:

    Coordinates_Distances: files with the distances between sampling sites and the coordinates of the sampling sites. 
        "OverWater_Comb12_Distance.csv" matrix of pairwise over-water distance (measured in Google Earth as the distance between two sampling sites that does not cross land), with sites 1 and 2 combined due to only one particle release site corresponding to both sites. Used in "Potential_Connectivity.R" script
     "OverWater_Distance.csv" matrix of pairwise over-water distance (measured in Google Earth as the distance between two sampling sites that does not cross land), used in "PRBI_IBD_Analysis.R" and "Potential_Connectivity.R" scripts
     "PRBI_AvgLatLong.csv" and "PRBI_Comb12_AvgLatLong" gives average latitude and longitude for individuals found at each sampling site, "PRBI_Comb12_AvgLatLong" combines populations 1 and 2 as run in the potential connectivity analyses.
     
     Fst_Matrices: Copies of Genepop pairwise Fst outputs in a .csv matrix format that can be loaded in RStudio. "2009-11-12 All Pop All Loci" pairwise Fst of all populations at all 16 loci. "2009-11-12_comb12_FstMatrix" pairwise Fst with populations 1 and 2 combined for potential connectivity analyses. 
     
   Genepop_Outputs: Outputs of Genepop pairwise Fst and Hardy-Weinberg functions. See "PRBI_IBD_Analysis.R" script for corresponding input files and code.
   
   Ne_Estimator: Inputs and output files from running NeEstimator (Do et al. 2014) for all sites in the IBD study region grouped together as 1 cohort with a pcrit of 0.02 and a monogamous mating model. "PRBI_genpop_1cohort_PopsExcluded.gen.txt" is the input file. "PRBI_genpop_1cohort_PopsExcluded.genNe.txt" gives the output Nb (Number of breeders) estimates and "PRBI_genpop_1cohort_PopsExcluded.genNoDat.txt" gives a summary of any loci with no data for any individuals.
   
   "PRBI_genepop_2009-11-12_comb12.gen.txt" genotypes for all Premnas biaculeatus individuals at 16 microsatellite loci in Genepop format. All individuals previously in sites 1 and 2 are combined into one population. Used as the input file for Fst calculations for the potential connectivity analyses.
     
R_Code: R scripts used for data analysis under the RProject file "Premnas-biaculeatus.RProj"

     "PRBI_IBD_Analysis.R" calculations of pairwise Fst matrices, testing for Hardy-Weinberg, running Mantel tests to test for IBD patterns, plots of geographic vs. genetic distance
     "Potential_Connectivity.R" identifies release sites for each sampling site from Thompson et al. 2018, tests for correlations between potential connectivity and genetic differentiation
     "Census_Density_Estimates.R" calculates census density for the IBD study region and bootstraps estimates to generate 95% confidence intervals
     "Estimate_Sigma" error propagation for estimates of density and IBD slope (m). Sources scripts from "Scripts for Sigma Estimates" folder.
     "other_study_dispersal_spreads" converts dispersal distances published for other coral reef species to dispersal spread for comparison purposes. Reported in the "Comparison to other coral reef species" section of the Discussion.
