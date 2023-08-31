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
  write_tsv('stickleback_maf0.05_ldpruned_filtered.tped', 
            col_names = F)
# 
snps = read_table('stickleback_maf0.05_ldpruned_filtered.tped', 
                  col_names = F)


df = snps %>% 
  dplyr::select(1, 
                5:length(snps)) %>% 
  filter(X1 == 'chr_XXI') %>% 
  dplyr::select(-X1) %>% 
  dplyr::select(-X5)

eigen = eigen_windows(df, 
                      win = 50, 
                      k = 5)
windist = pc_dist(eigen, 
                  npc = 5) %>% 
  as_tibble()

window_data = snps %>% 
  select(1:4) %>% 
  filter(X1 == 'chr_XXI') %>% 
  mutate(window = ceiling(row_number()/50)) %>% 
  group_by(window) %>% 
  mutate(mean_window = mean(X4)) %>% 
  distinct(mean_window, 
           .keep_all = T) %>% 
  filter(window %in% 1:nrow(windist))

combo_data = bind_cols(window_data, 
                       windist)


MDS_data = cmdscale(combo_data[7:length(combo_data)], 
                    eig = TRUE, 
                    k = 5)

MDS_points = MDS_data$points %>% 
  as_tibble() %>% 
  rename(MDS_Points1 = V1, 
         MDS_Points2 = V2, 
         MDS_points3 = V3, 
         MDS_points4 = V4, 
         MDS_points5 = V5)

combo_data = bind_cols(combo_data, 
                       MDS_points) 

