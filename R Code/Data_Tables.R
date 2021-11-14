##Data Tables for Manuscript##

library(gt)

#Sample size table

sample_size <- read_csv("Sample_Sizes.csv")

sample_size <- sapply(sample_size, as.character)          # Convert all columns to character

sample_size[is.na(sample_size)] <- ""   # Replace NA with blank

sample_size <- sapply(sample_size, as.numeric)

sample_size %>%
  gt() %>%
  tab_header(title="Table 2. Sample sizes of P. biaculeatus at each site.")

