##########################################
## Divergence in immune genes
##
## 04.11.2025
##
## Matthew K. Brachmann (phdmattyb)
#########################################

setwd('~/Parsons_Postdoc/Stickleback_Genomic/vcf_filter/')


library(tidyverse)


ASHN_Fst = read_tsv('ASHN_Fst_values.fst') %>% 
  na.omit() %>%  ##pull out na's
  mutate(FST_zero = if_else(FST < 0, 0, FST))
MYV_Fst = read_tsv('MYV_Fst_values.fst')%>% 
  na.omit() %>%  ##pull out na's
  mutate(FST_zero = if_else(FST < 0, 0, FST))
SKR_Fst = read_tsv('SKR_Fst_values.fst')%>% 
  na.omit() %>%  ##pull out na's
  mutate(FST_zero = if_else(FST < 0, 0, FST))
GTS_CSWY_Fst = read_tsv('GTS_CSWY_Fst_values.fst')%>% 
  na.omit() %>%  ##pull out na's
  mutate(FST_zero = if_else(FST < 0, 0, FST))
WC_Fst = read_tsv('Warm_Cold_Fst.fst') %>% 
  na.omit() %>% 
  mutate(FST_zero = if_else(FST < 0, 0, FST))



# spi1b gene --------------------------------------------------------------

#chrII:9,234,445..9,240,960

ASHN_Fst %>% 
  filter(CHR == 'chr_II') %>% 
  filter(POS >= 9234445, 
         POS <= 9240960)
MYV_Fst %>% 
  filter(CHR == 'chr_II') %>% 
  filter(POS >= 9234445, 
         POS <= 9240960)
SKR_Fst %>% 
  filter(CHR == 'chr_II') %>% 
  filter(POS >= 9234445, 
         POS <= 9240960)
GTS_CSWY_Fst %>% 
  filter(CHR == 'chr_II') %>% 
  filter(POS >= 9234445, 
         POS <= 9240960)
WC_Fst %>% 
  filter(CHR == 'chr_II') %>% 
  filter(POS >= 9234445, 
         POS <= 9240960)

