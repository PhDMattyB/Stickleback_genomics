##############################
## Per pop outlier gene names
##
## Matt Brachmann (PhDMattyB)
##
## 14.04.2026
##
##############################

setwd("~/Parsons_Postdoc/Stickleback_Genomic/Stickleback_Annotation_features/")

library(tidyverse)


ASHN_gene_names = read_tsv('ASHN_FST_0.5%_outlier_gene_names_only_NEW.tsv', 
                           col_names = 'GeneID') %>% 
  mutate(Population = 'ASHN')


MYV_gene_names = read_tsv('MYV_FST_0.5%_outlier_gene_names_only_NEW.tsv', 
                           col_names = 'GeneID') %>% 
  mutate(Population = 'MYV')

SKR_gene_names = read_tsv('SKR_FST_0.5%_outlier_gene_names_only_NEW.tsv', 
                           col_names = 'GeneID') %>% 
  mutate(Population = 'SKR')

GTS_CSWY_gene_names = read_tsv('GTS_CSWY_FST_0.5%_outlier_gene_names_only_NEW.tsv', 
                           col_names = 'GeneID') %>% 
  mutate(Population = 'GTS-GAR')


gene_name_list = bind_rows(ASHN_gene_names, 
                           MYV_gene_names, 
                           SKR_gene_names, 
                           GTS_CSWY_gene_names)

gene_name_list %>% 
  write_csv('NatComms_SuppFigureS2-3_data.csv')
