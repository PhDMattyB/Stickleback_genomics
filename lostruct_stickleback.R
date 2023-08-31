##############################
## lostruct stickleback chr21
##
## Matt Brachmann (PhDMattyB)
##
## 30.08.2023
##
##############################

setwd('~/Parsons_Postdoc/Stickleback_Genomic/vcf_filter/')

library(tidyverse)
library(data.table)

ped_data = read_table('stickleback_maf0.05_ldpruned_filtered.ped',
                      col_names = c('PopulationID', 
                                    'IndividualID', 
                                    'MaternalID', 
                                    'PaternalID', 
                                    'Sex', 
                                    'Phenotype', 
                                    map_data$SNP))

map_data = read_tsv('stickleback_maf0.05_ldpruned_filtered.map', 
                    col_names = c('CHR', 
                                  'SNP', 
                                  'GPOS', 
                                  'POS'))

Create_tped(ped = ped_data, 
            map = map_data)
