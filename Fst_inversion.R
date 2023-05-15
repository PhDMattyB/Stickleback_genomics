##############################
## FST inversion hunting
##
## Matt Brachmann (PhDMattyB)
##
## 03.05.2023
##
##############################


setwd('~/Parsons_Postdoc/Stickleback_Genomic/vcf_filter/')


library(tidyverse)


# Fst set up Plink ---------------------------------------------------------------

ped_test = read_table2('stickleback_maf0.05_ldpruned_filtered.ped', 
                       col_names = F)

identifiers = read_csv('stickleback_identifiers.csv')
identifiers = mutate(.data = identifiers,
                     type = as.factor(case_when(
                       population == 'ASHNC' ~ 'Cold',
                       population == 'ASHNW' ~ 'Warm',
                       population == 'CSWY' ~ 'Manmade',
                       population == 'GTS' ~ 'Thermal',
                       population == 'MYVC' ~ 'Cold',
                       population == 'MYVW' ~ 'Warm',
                       population == 'SKRC' ~ 'Cold',
                       population == 'SKRW' ~ 'Warm')))

ped_ids = read_table2('stickleback_maf0.05_ldpruned_filtered.fam', 
                      col_names = F) %>%
  dplyr::select(X1,
                X2)
ped_ids = bind_cols(ped_ids, 
                    identifiers)


## Need to split based on GTS vs MYVC
## Need to split based on GTS vs SKRC, 
## Need to split based on GTS cs ASHNC

## Need to split based on CSWY vs MYVW
## Need to split based on CSWY vs SKRW
## Need to split based on CSWY vs ASHNW

ped_ids %>%
  filter(population %in% c('CSWY',
                         'ASHNW')) %>%
  select(X1,
         X2,
         type) %>%
  rename(`#population` = 1,
         individual_id = 2) %>%
  write_tsv('CSWY_ASHNW_Fst_grouping.txt')


## Need to make a ped and map file for each of these comparisons
## the ped file is waaaay to big to open in R
## use the --keep or --keep-fam flags in plink to filter the 
## populations out. 
## the --keep file needs to be a text file with family and individual
## identifiers


ped_ids %>% 
  filter(population %in% c('CSWY', 
                           'ASHNW')) %>% 
  # filter(population %in% c('GTS', 
  #                          'CSWY')) %>% 
  # filter(type %in% c('Warm', 
  #                          'Cold')) %>%
  dplyr::select(X1, 
                X2) %>% 
  # rename(`#population` = population, 
  # individual_ID = X1) %>% 
  write_tsv('CSWY_ASHNW_keep.txt', 
            col_names = F)

##


# GTS vs cold pops FST ---------------------------------------------------------------

GTS_MYVC = read_tsv('GTS_MYVC_FST.fst') %>% 
  na.omit() %>%  ##pull out na's
  mutate(FST_zero = if_else(FST < 0, 0, FST))%>% 
  stickle_CHR_reorder() %>% 
  dist_cal()

GTS_SKRC = read_tsv('GTS_SKRC_FST.fst') %>% 
  na.omit() %>%  ##pull out na's
  mutate(FST_zero = if_else(FST < 0, 0, FST))%>% 
  stickle_CHR_reorder() %>% 
  dist_cal()

GTS_ASHNC = read_tsv('GTS_ASHNC_FST.fst') %>% 
  na.omit() %>%  ##pull out na's
  mutate(FST_zero = if_else(FST < 0, 0, FST))%>% 
  stickle_CHR_reorder() %>% 
  dist_cal()


