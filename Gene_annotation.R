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


read_tsv('stickleback_v5_ensembl_genes.gff3.gz', 
         col_names = F, 
         skip = 8) %>% 
  select(X3) %>% 
  distinct()

genes_start = read_tsv('stickleback_v5_ensembl_genes.gff3.gz', 
                  col_names = F, 
                  skip = 8) %>% 
  # filter(X3 %in% c('gene', 
  #                  'exon', 
  #                  'CDS')) %>% 
  group_by(X1) %>% 
  arrange(X4, 
          X5) %>% 
  ## arrange each gene by its start and end points on each chromosome
  mutate(mid = X4 + (X5-X4)/2) %>% 
  dplyr::select(X1, 
                X3:X5, 
                X9:mid) %>% 
  rename(chromosome = X1, 
         feature = X3, 
         start = X4, 
         end = X5, 
         gene_id = X9, 
         position = mid) %>% 
  dplyr::select(chromosome, 
                feature, 
                start, 
                gene_id) %>% 
  rename(position = start)

genes_end = read_tsv('stickleback_v5_ensembl_genes.gff3.gz', 
                       col_names = F, 
                       skip = 8) %>% 
  # filter(X3 %in% c('gene', 
  #                  'exon', 
  #                  'CDS')) %>% 
  group_by(X1) %>% 
  arrange(X4, 
          X5) %>% 
  ## arrange each gene by its start and end points on each chromosome
  mutate(mid = X4 + (X5-X4)/2) %>% 
  dplyr::select(X1, 
                X3:X5, 
                X9:mid) %>% 
  rename(chromosome = X1, 
         feature = X3, 
         start = X4, 
         end = X5, 
         gene_id = X9, 
         position = mid) %>% 
  dplyr::select(chromosome, 
                feature, 
                end, 
                gene_id) %>% 
  rename(position = end)

# %>% 
#   ## calculate the mid point of each gene from the start and end points
#   # filter(mid >= 406599, 
#   #        mid <= 33118168) %>%
#   pull(X9) %>% 
#   str_split_fixed(';', 
#                   n = 5) %>% 
#   as_tibble() %>%
#   dplyr::select(V3) %>% 
#   separate(col = V3, 
#            into = c('geneid', 
#                     'accenssionid'), 
#            sep = ',') %>% 
#   na.omit() %>% 
#   dplyr::select(geneid) %>% 
#   separate(col = geneid, 
#            into = c('garbage', 
#                     'geneid'), 
#            sep = ':') %>% 
#   dplyr::select(geneid)

methy_outliers = read_tsv('~/Parsons_Postdoc/Stickleback_Genomic/Methylation_outliers.txt')%>% 
  separate(col = TETWarm, 
           into = c('chromosome', 
                    'position'),
           sep = '-') %>% 
  group_by(chromosome) %>% 
  arrange(position)

methy_outliers$position = as.numeric(methy_outliers$position)


genes_start_lineup = inner_join(methy_outliers, 
           genes_start, 
           by = c('chromosome', 
                  'position')) %>% 
  View()

genes_end_lineup = inner_join(methy_outliers, 
           genes_end, 
           by = c('chromosome', 
                  'position')) 

genes_end_lineup %>% 
  write_tsv('Methylation_data_genes_lineup.txt')
  

# Alternative assembly information ----------------------------------------

read_tsv('genomic.gff', 
         col_names = F, 
         skip = 8) %>% 
  select(X3) %>% 
  distinct()


genes_start = read_tsv('genomic.gff', 
                       col_names = F, 
                       skip = 8) %>% 
  # filter(X3 %in% c('gene', 
  #                  'exon', 
  #                  'CDS')) %>% 
  group_by(X1) %>% 
  arrange(X4, 
          X5) %>% 
  ## arrange each gene by its start and end points on each chromosome
  mutate(mid = X4 + (X5-X4)/2) %>% 
  dplyr::select(X1, 
                X3:X5, 
                X9:mid) %>% 
  rename(chromosome = X1, 
         feature = X3, 
         start = X4, 
         end = X5, 
         gene_id = X9, 
         position = mid) %>% 
  dplyr::select(chromosome, 
                feature, 
                start, 
                gene_id) %>% 
  rename(position = start)

genes_start %>% 
  dplyr::select(chromosome) %>% 
  distinct() %>%
  arrange(chromosome) %>% 
  View()
## there is a weird chromosome thats not in the genome annotation details
## Also there is a weird set of chromosomes with a hashtag

methy_outliers = read_tsv('~/Parsons_Postdoc/Stickleback_Genomic/Methylation_outliers.txt')%>% 
  separate(col = TETWarm, 
           into = c('chromosome', 
                    'position'),
           sep = '-') %>% 
  group_by(chromosome) %>% 
  arrange(position)

methy_outliers$position = as.numeric(methy_outliers$position)


genes_start_lineup = inner_join(methy_outliers, 
                                genes_start, 
                                by = c('chromosome', 
                                       'position')) %>% 
  View()

genes_end_lineup = inner_join(methy_outliers, 
                              genes_end, 
                              by = c('chromosome', 
                                     'position')) 

genes_end_lineup %>% 
  write_tsv('Methylation_data_genes_lineup.txt')

