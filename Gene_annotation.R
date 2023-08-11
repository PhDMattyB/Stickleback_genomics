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


gene_annotation = read_tsv('stickleback_v5_ensembl_genes.gff3.gz', 
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
  na.omit()

methy_outliers = read_tsv('~/Parsons_Postdoc/Stickleback_Genomic/Methylation_outliers.txt')%>% 
  separate(col = TETWarm, 
           into = c('chromosome', 
                    'position'),
           sep = '-') %>% 
  group_by(chromosome) %>% 
  arrange(position)

methy_outliers$position = as.numeric(methy_outliers$position)

## create a 1Kb window around the snp identified
## settting the start and end range for the methylation data
## this number can be as big or small as you want depending
## on how liberal you want it

methy_out_window = methy_outliers %>%
  group_by(chromosome) %>% 
  mutate(start = position-100, 
         end = position+100)

library(data.table)
  
setDT(methy_out_window)
setDT(gene_annotation)

setkey(methy_out_window, 
       chromosome, 
       start, 
       end)

gene_overlap = foverlaps(gene_annotation, 
          methy_out_window, 
          type="any")


gene_overlap_tib = as_tibble(gene_overlap) %>% 
  na.omit() 


gene_name_1 = gene_overlap_tib %>% 
  # pull(gene_id) %>% 
  as_tibble() %>% 
  dplyr::select(chromosome, 
                position, 
                start, 
                end, 
                gene_id, 
                F118, 
                F218, 
                TETWarm.F118, 
                TETWarm.F218, 
                TETWarm.F118.F218) %>% 
  separate(col = gene_id, 
           into = c('ensemble_id', 
                    'gene_name', 
                    'parent_code', 
                    'gene_name2'), 
           sep = ';') %>%
  separate(col = gene_name, 
           into = c('Garbage', 
                    'gene_name'), 
           sep = '=') %>% 
  dplyr::select(chromosome, 
                position, 
                start, 
                end, 
                gene_name) %>% 
  na.omit()

gene_name_2 = gene_overlap_tib %>% 
  # pull(gene_id) %>% 
  as_tibble() %>% 
  dplyr::select(chromosome, 
                position, 
                start, 
                end, 
                gene_id,
                F118, 
                F218, 
                TETWarm.F118, 
                TETWarm.F218, 
                TETWarm.F118.F218) %>% 
  separate(col = gene_id, 
           into = c('ensemble_id', 
                    'gene_name', 
                    'parent_code', 
                    'gene_name2'), 
           sep = ';') %>%
  separate(col = parent_code, 
           into = c('Garbage', 
                    'gene_name'), 
           sep = '=') %>% 
  dplyr::select(chromosome, 
                position, 
                start, 
                end, 
                gene_name) %>% 
  na.omit()


methy_genes = bind_rows(gene_name_1, 
                        gene_name_2) %>% 
  arrange(chromosome, 
          position) %>% 
  distinct(gene_name, 
           .keep_all = T) %>% 
  filter(!grepl('ENSG', 
                gene_name))
  
methy_genes %>% 
  write_csv('Methylation_outlier_genes.csv')
