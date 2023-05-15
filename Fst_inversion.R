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

axisdf_MYVC = axis_df(GTS_MYVC)
axisdf_SKRC = axis_df(GTS_SKRC)
axisdf_ASHNC = axis_df(GTS_ASHNC)

GTS_MYVC_outs = GTS_MYVC[GTS_MYVC$FST_zero > quantile(GTS_MYVC$FST_zero, 
                                  prob = 1-5/100),]

GTS_SKRC_outs = GTS_SKRC[GTS_SKRC$FST_zero > quantile(GTS_SKRC$FST_zero, 
                                                      prob = 1-5/100),]


GTS_ASHNC_outs = GTS_ASHNC[GTS_ASHNC$FST_zero > quantile(GTS_ASHNC$FST_zero, 
                                                      prob = 1-5/100),]

Fst_manhattan(non_outs = GTS_MYVC, 
              outs = GTS_MYVC_outs,
              axisdf = axisdf_MYVC, 
              xval = BPcum, 
              yval = FST_zero, 
              chr = GTS_MYVC$CHR, 
              out_col = '#d62828',)

Fst_manhattan(non_outs = GTS_ASHNC, 
              outs = GTS_ASHNC_outs, 
              axisdf = axisdf_ASHNC, 
              xval = BPcum, 
              yval = FST_zero, 
              chr = GTS_ASHNC$CHR, 
              out_col = '#06d6a0')

Fst_manhattan(non_outs = GTS_SKRC, 
              outs = GTS_SKRC_outs, 
              axisdf = axisdf_SKRC, 
              xval = BPcum, 
              yval = FST_zero, 
              chr = GTS_SKRC$CHR, 
              out_col = '#5f0f40')
