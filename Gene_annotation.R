##############################
## Stickleback gene annotation 
##
## Matt Brachmann (PhDMattyB)
##
## 08.08.2023
##
##############################


setwd('~/Parsons_Postdoc/Stickleback_Genomic/Stickleback_Annotation_features/')

library(tidyverse)


read_tsv('genomic.gff', 
         col_names = F, 
         skip = 8) %>% 
  select(X3) %>% 
  distinct()

geneid = read_tsv('genomic.gff', 
                  col_names = F, 
                  skip = 8) %>% 
  filter(X3 %in% c('gene', 
                   'exon', 
                   'CDS')) %>% 
  group_by(X1) %>% 
  arrange(X4, 
          X5) %>% 
  ## arrange each gene by its start and end points on each chromosome
  mutate(mid = X4 + (X5-X4)/2) %>% 
  ## calculate the mid point of each gene from the start and end points
  # filter(mid >= 406599, 
  #        mid <= 33118168) %>%
  pull(X9) %>% 
  str_split_fixed(';', 
                  n = 5) %>% 
  as_tibble() %>%
  dplyr::select(V3) %>% 
  separate(col = V3, 
           into = c('geneid', 
                    'accenssionid'), 
           sep = ',') %>% 
  na.omit() %>% 
  dplyr::select(geneid) %>% 
  separate(col = geneid, 
           into = c('garbage', 
                    'geneid'), 
           sep = ':') %>% 
  dplyr::select(geneid)
