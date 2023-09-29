##############################
## Grandparent genotypes for Corin
##
## Matt Brachmann (PhDMattyB)
##
## 29.09.2023
##
##############################


setwd('~/Parsons_Postdoc/Stickleback_Genomic/vcf_filter/SKR_corrected/')

library(tidyverse)
library(data.table)

map_data = read_tsv('stickleback_maf0.05_ldpruned_filtered.map', 
                    col_names = c('CHR', 
                                  'SNP', 
                                  'GPOS', 
                                  'POS'))

ped_data = read_table('stickleback_maf0.05_ldpruned_filtered.ped',
                      col_names = c('PopulationID', 
                                    'IndividualID', 
                                    'MaternalID', 
                                    'PaternalID', 
                                    'Sex', 
                                    'Phenotype', 
                                    map_data$SNP))


Corin_skr_map = read_tsv('SKR_plink.map', 
                         col_names = c('CHR', 
                                       'SNP', 
                                       'GPOS', 
                                       'POS')) %>% 
  select(-SNP)

## SNP ids are messed up, need to make them myself

sep_map = Corin_skr_map %>% 
  dplyr::select(CHR, 
                POS) %>% 
  separate(col = CHR, 
           into = c('label', 
                    'chr_num'), 
           sep = "r")
chr_label = rep('chr', 
                nrow(sep_map)) %>% 
  as_tibble()

snp_ID = bind_cols(sep_map, 
          chr_label) %>% 
  select(value,
         chr_num, 
         POS) %>% 
  unite(col = 'SNP', 
        c('value', 
          'chr_num', 
          'POS'), 
        sep = '_')

Corin_skr_map = bind_cols(Corin_skr_map, 
          snp_ID) %>% 
  select(CHR, 
         SNP, 
         GPOS, 
         POS)

# Corin_skr_map %>%
#   group_by(CHR) %>%
#   summarize(n = n()) %>%
#   View()

Corin_skr_ped = read_table('SKR_plink.ped',
                           col_names = c('PopulationID', 
                                         'IndividualID', 
                                         'MaternalID', 
                                         'PaternalID', 
                                         'Sex', 
                                         'Phenotype', 
                                         Corin_skr_map$SNP))
Corin_skr_ped %>% 
  select(PopulationID) %>% 
  View()

Corin_test_snp = Corin_skr_ped %>% 
  select(starts_with('chr_'))%>%
  mutate(across(everything(), 
                as.character))

whole_genome_snp_test = ped_data %>%
  select(PopulationID, 
         starts_with('chr_')) %>%
  filter(PopulationID %in% c('Sample_54-SKRW_4', 
                             'Sample_38-SKRC_4')) %>% 
  select(-PopulationID) %>% 
  mutate(across(everything(), 
                as.character))

## fuck the inner join didn't work
overlapping_snps = inner_join(whole_genome_snp_test, 
           Corin_test_snp)

Corin_skr_map 

map_data

inner_join(map_data %>% 
           Corin_skr_map, 
           by = 'SNP')

