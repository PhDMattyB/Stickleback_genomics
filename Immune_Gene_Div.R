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

ASHN_spi1b = ASHN_Fst %>% 
  filter(CHR == 'chr_II') %>% 
  filter(POS >= 9234445, 
         POS <= 9240960) %>% 
  mutate(Population = 'ASHN', 
         Gene = 'spi1b')
MYV_spi1b = MYV_Fst %>% 
  filter(CHR == 'chr_II') %>% 
  filter(POS >= 9234445, 
         POS <= 9240960) %>% 
  mutate(Population = 'MYV', 
         Gene = 'spi1b')
SKR_spi1b = SKR_Fst %>% 
  filter(CHR == 'chr_II') %>% 
  filter(POS >= 9234445, 
         POS <= 9240960) %>% 
  mutate(Population = 'SKR', 
         Gene = 'spi1b')
GTS_CSWY_spi1b = GTS_CSWY_Fst %>% 
  filter(CHR == 'chr_II') %>% 
  filter(POS >= 9234445, 
         POS <= 9240960) %>% 
  mutate(Population = 'GTS-GAR', 
         Gene = 'spi1b')
WC_spi1b = WC_Fst %>% 
  filter(CHR == 'chr_II') %>% 
  filter(POS >= 9234445, 
         POS <= 9240960) %>% 
  mutate(Population = 'All Warm-Cold', 
         Gene = 'spi1b')

spi1b = bind_rows(ASHN_spi1b, 
          MYV_spi1b, 
          SKR_spi1b, 
          GTS_CSWY_spi1b, 
          WC_spi1b)


# stat6 -------------------------------------------------------------------
#chrXII:10,262,569..10,274,676

ASHN_stat6 = ASHN_Fst %>% 
  filter(CHR == 'chr_XII') %>% 
  filter(POS >= 10262569, 
         POS <= 10274678) %>% 
  mutate(Population = 'ASHN', 
         Gene = 'stat6')
MYV_stat6 = MYV_Fst %>% 
  filter(CHR == 'chr_XII') %>% 
  filter(POS >= 10262569, 
         POS <= 10274678) %>% 
  mutate(Population = 'MYV', 
         Gene = 'stat6')
SKR_stat6 = SKR_Fst %>% 
  filter(CHR == 'chr_XII') %>% 
  filter(POS >= 10262569, 
         POS <= 10274678) %>% 
  mutate(Population = 'SKR', 
         Gene = 'stat6')
GTS_CSWY_stat6 = GTS_CSWY_Fst %>% 
  filter(CHR == 'chr_XII') %>% 
  filter(POS >= 10262569, 
         POS <= 10274678) %>% 
  mutate(Population = 'GTS-GAR', 
         Gene = 'stat6')
WC_stat6 = WC_Fst %>% 
  filter(CHR == 'chr_XII') %>% 
  filter(POS >= 10262569, 
         POS <= 10274678) %>% 
mutate(Population = 'All Warm-Cold', 
       Gene = 'stat6')


stat6 = bind_rows(ASHN_stat6, 
          MYV_stat6, 
          SKR_stat6, 
          GTS_CSWY_stat6, 
          WC_stat6)

# cyp3a48 -----------------------------------------------------------------
#chrXII:12,613,328..12,615,383
#chrXII:8,154,600..8,158,564 


ASHN_cyp3a48 = ASHN_Fst %>% 
  filter(CHR == 'chr_XII') %>% 
  filter(POS >= 8154600, 
         POS <= 8158564)
MYV_cyp3a48 = MYV_Fst%>% 
  filter(CHR == 'chr_XII') %>% 
  filter(POS >= 8154600, 
         POS <= 8158564)
SKR_cyp3a48 = SKR_Fst %>% 
  filter(CHR == 'chr_XII') %>% 
  filter(POS >= 8154600, 
         POS <= 8158564)
GTS_CSWY_cyp3a48 = GTS_CSWY_Fst%>% 
  filter(CHR == 'chr_XII') %>% 
  filter(POS >= 8154600, 
         POS <= 8158564)
WC_cyp3a48 = WC_Fst %>% 
  filter(CHR == 'chr_XII') %>% 
  filter(POS >= 8154600, 
         POS <= 8158564)

cyp3a48 = bind_rows(ASHN_cyp3a48, 
                    MYV_cyp3a48, 
                    SKR_cyp3a48, 
                    GTS_CSWY_cyp3a48, 
                    WC_cyp3a48)

# ptpn6 -------------------------------------------------------------------

##chrXX:7,807,984..7,819,960

ASHN_ptpn6 = ASHN_Fst %>% 
  filter(CHR == 'chr_XX') %>% 
  filter(POS >= 7807984, 
         POS <= 7819960)%>% 
  mutate(Population = 'ASHN', 
         Gene = 'ptpn6')
MYV_ptpn6 = MYV_Fst %>% 
  filter(CHR == 'chr_XX') %>% 
  filter(POS >= 7807984, 
         POS <= 7819960)%>% 
  mutate(Population = 'MYV', 
         Gene = 'ptpn6')
SKR_ptpn6 = SKR_Fst %>% 
  filter(CHR == 'chr_XX') %>% 
  filter(POS >= 7807984, 
         POS <= 7819960) %>% 
  mutate(Population = 'SKR', 
         Gene = 'ptpn6')
GTS_CSWY_ptpn6 = GTS_CSWY_Fst %>% 
  filter(CHR == 'chr_XX') %>% 
  filter(POS >= 7807984, 
         POS <= 7819960)%>% 
  mutate(Population = 'GTS-GAR', 
         Gene = 'ptpn6')
WC_ptpn6 = WC_Fst %>% 
  filter(CHR == 'chr_XX') %>% 
  filter(POS >= 7807984, 
         POS <= 7819960)%>% 
  mutate(Population = 'All Warm-Cold', 
         Gene = 'ptpn6')

ptpn6 = bind_rows(ASHN_ptpn6, 
                  MYV_ptpn6, 
                  SKR_ptpn6, 
                  GTS_CSWY_ptpn6, 
                  WC_ptpn6)


# all gene info ----------------------------------------------------------

bind_rows(spi1b, 
          stat6,
          cyp3a48,
          ptpn6) %>% 
  write_tsv('Immune_Gene_FST_info.txt')
