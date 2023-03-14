##############################
## 
## selective sweep scan
##
## Matt Brachmann (PhDMattyB)
##
## 14.03.2023
##
##############################

setwd('~/Parsons_Postdoc/Stickleback_Genomic/vcf_filter/')
library(tidyverse)
library(devtools)

# zalpha package ----------------------------------------------------------

# install_github("chorscroft/zalpha")

library(zalpha)

## Need a vector of snp locations
stickle_map = read_tsv('stickleback_maf0.05_ldpruned_filtered.map', 
                       col_names = c('CHR', 
                                     'SNP', 
                                     'GENO_POS', 
                                     'POS'))
snp_loc = stickle_map %>% 
  select(POS)

## This is really not working and R keep aborting
# stickle_ped = read_table2('stickleback_maf0.05_ldpruned_filtered.ped') %>% 
#   dplyr::select(1:13)
