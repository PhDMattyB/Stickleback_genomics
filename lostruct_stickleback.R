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

ped_data = read_table('stickleback_maf0.05_ldprunded_filtered.raw',
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
# map_data %>%
#   # select(CHR) %>% 
#   # distinct() %>% 
#   # View()
#   filter(CHR == 'chr_XXI')


tped = Create_tped(ped = ped_data, 
            map = map_data)

tped %>% 
  write_tsv('stickleback_maf0.05_ldpruned_filtered.tped')

snps = read_tped('stickleback_maf0.05_ldpruned_filtered.tped')


df = tped %>% 
  dplyr::select(CHR, 
                5:length(tped)) %>% 
  filter(CHR == 'chr_XXI') %>% 
  dplyr::select(-CHR) %>% 
  dplyr::select(-IID) %>% 
  as.data.frame()

eigen = eigen_windows(df, 
                      win = 50, 
                      k = 5)

lostruct_run(data = tped, 
             chr = 'chr_XXI', 
             window_size = 50, 
             k_value = 5)

eigen_windows()
