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
library(data.table)

# 
# read_tsv('stickleback_v5_ensembl_genes.gff3.gz', 
#          col_names = F, 
#          skip = 1) %>% 
#   select(X3) %>% 
#   distinct()


gene_annotation = read_tsv('stickleback_v5_ensembl_genes.gff3.gz', 
                     col_names = F, 
                     skip = 1) %>% 
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

gene_annotation %>% 
  arrange(chromosome) %>%
  # filter(feature == 'gene')
  group_by(feature) %>% 
  distinct(feature)




#### ASHN outliers ----------------------------------------------------
  
ASHN_div_snps = read_csv('ASHN_TOP_DAWG_Fst_clean.csv') %>% 
  select(-NMISS, 
         -FST) %>%
  rename(position = POS, 
         chromosome = CHR, 
         FST = FST_zero) %>% 
  filter(value == 'Outlier') %>% 
  arrange(chromosome, 
          position) %>% 
  group_by(chromosome)
## create a 100bp window around the snp identified
## settting the start and end range for the methylation data
## this number can be as big or small as you want depending
## on how liberal you want it

ASHN_div_window = ASHN_div_snps %>%
  group_by(chromosome)%>%
  mutate(start = position-100,
         end = position+100) %>% 
  separate(col = chromosome, 
           into = c('chr', 
                    'chr_name'), 
           sep = '_') %>% 
  unite(chromosome, 
        chr:chr_name,
        sep = '',
        remove = F) %>% 
  select(-chr, 
         -chr_name, 
         -SNP) %>% 
  group_by(chromosome)


setDT(ASHN_div_window)
setDT(gene_annotation)

setkey(ASHN_div_window, 
       chromosome, 
       start, 
       end)

ASHN_gene_overlap = foverlaps(gene_annotation,
                         ASHN_div_window,
                         # by.x = start,
                         # by.y = end,
                         type="any")


ASHN_gene_overlap_tib = as_tibble(ASHN_gene_overlap) %>% 
  na.omit() %>% 
  filter(chromosome != 'chrUn') %>% 
  arrange(chromosome, 
          position)

gene_name_1 = ASHN_gene_overlap_tib %>% 
  # pull(gene_id) %>% 
  as_tibble() %>% 
  dplyr::select(chromosome, 
                position, 
                start,
                i.start, 
                end, 
                i.end,
                gene_id,
                feature,
                FST) %>% 
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
                i.start,
                end, 
                i.end,
                FST,
                feature,
                gene_name) %>% 
  na.omit()

gene_name_2 = ASHN_gene_overlap_tib %>% 
  # pull(gene_id) %>% 
  as_tibble() %>% 
  dplyr::select(chromosome, 
                position, 
                start,
                i.start,
                end, 
                i.end,
                gene_id,
                feature,
                FST,) %>% 
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
                i.start,
                end, 
                i.end,
                FST,
                feature,
                gene_name) %>% 
  na.omit()


ASHN_FST_out_genes = bind_rows(gene_name_1, 
                        gene_name_2) %>% 
  arrange(chromosome, 
          position) %>% 
  distinct(gene_name, 
           .keep_all = T) %>% 
  filter(!grepl('ENSG', 
                gene_name)) %>% 
  filter(feature == 'gene')

ASHN_FST_out_genes %>% 
  write_csv('ASHN_FST_0.5%_outlier_genes.csv')

ASHN_FST_out_genes %>% 
  select(gene_name) %>% 
  write_tsv('ASHN_FST_0.5%_outlier_gene_names_only.tsv', 
            col_names = F)


#### MYV outliers ----------------------------------------------------

MYV_div_snps = read_csv('MYV_TOP_DAWG_Fst_clean.csv') %>% 
  select(-NMISS, 
         -FST) %>%
  rename(position = POS, 
         chromosome = CHR, 
         FST = FST_zero) %>% 
  filter(value == 'Outlier') %>% 
  arrange(chromosome, 
          position) %>% 
  group_by(chromosome)
## create a 100bp window around the snp identified
## settting the start and end range for the methylation data
## this number can be as big or small as you want depending
## on how liberal you want it

MYV_div_window = MYV_div_snps %>%
  group_by(chromosome)%>%
  mutate(start = position-100,
         end = position+100) %>% 
  separate(col = chromosome, 
           into = c('chr', 
                    'chr_name'), 
           sep = '_') %>% 
  unite(chromosome, 
        chr:chr_name,
        sep = '',
        remove = F) %>% 
  select(-chr, 
         -chr_name, 
         -SNP) %>% 
  group_by(chromosome)


setDT(MYV_div_window)
setDT(gene_annotation)

setkey(MYV_div_window, 
       chromosome, 
       start, 
       end)

MYV_gene_overlap = foverlaps(gene_annotation,
                              MYV_div_window,
                              # by.x = start,
                              # by.y = end,
                              type="any")


MYV_gene_overlap_tib = as_tibble(MYV_gene_overlap) %>% 
  na.omit() %>% 
  filter(chromosome != 'chrUn') %>% 
  arrange(chromosome, 
          position)

gene_name_1 = MYV_gene_overlap_tib %>% 
  # pull(gene_id) %>% 
  as_tibble() %>% 
  dplyr::select(chromosome, 
                position, 
                start,
                i.start, 
                end, 
                i.end,
                gene_id,
                feature,
                FST) %>% 
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
                i.start,
                end, 
                i.end,
                FST,
                feature,
                gene_name) %>% 
  na.omit()

gene_name_2 = MYV_gene_overlap_tib %>% 
  # pull(gene_id) %>% 
  as_tibble() %>% 
  dplyr::select(chromosome, 
                position, 
                start,
                i.start,
                end, 
                i.end,
                gene_id,
                feature,
                FST,) %>% 
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
                i.start,
                end, 
                i.end,
                FST,
                feature,
                gene_name) %>% 
  na.omit()


MYV_FST_out_genes = bind_rows(gene_name_1, 
                               gene_name_2) %>% 
  arrange(chromosome, 
          position) %>% 
  distinct(gene_name, 
           .keep_all = T) %>% 
  filter(!grepl('ENSG', 
                gene_name)) %>% 
  filter(feature == 'gene')

MYV_FST_out_genes %>% 
  write_csv('MYV_FST_0.5%_outlier_genes.csv')

MYV_FST_out_genes %>% 
  select(gene_name) %>% 
  write_tsv('MYV_FST_0.5%_outlier_gene_names_only.tsv', 
            col_names = F)


#### SKR outliers ----------------------------------------------------

SKR_div_snps = read_csv('SKR_TOP_DAWG_Fst_clean.csv') %>% 
  select(-NMISS, 
         -FST) %>%
  rename(position = POS, 
         chromosome = CHR, 
         FST = FST_zero) %>% 
  filter(value == 'Outlier') %>% 
  arrange(chromosome, 
          position) %>% 
  group_by(chromosome)
## create a 100bp window around the snp identified
## settting the start and end range for the methylation data
## this number can be as big or small as you want depending
## on how liberal you want it

SKR_div_window = SKR_div_snps %>%
  group_by(chromosome)%>%
  mutate(start = position-100,
         end = position+100) %>% 
  separate(col = chromosome, 
           into = c('chr', 
                    'chr_name'), 
           sep = '_') %>% 
  unite(chromosome, 
        chr:chr_name,
        sep = '',
        remove = F) %>% 
  select(-chr, 
         -chr_name, 
         -SNP) %>% 
  group_by(chromosome)

setDT(SKR_div_window)
setDT(gene_annotation)

setkey(SKR_div_window, 
       chromosome, 
       start, 
       end)

SKR_gene_overlap = foverlaps(gene_annotation,
                              SKR_div_window,
                              # by.x = start,
                              # by.y = end,
                              type="any")


SKR_gene_overlap_tib = as_tibble(SKR_gene_overlap) %>% 
  na.omit() %>% 
  filter(chromosome != 'chrUn') %>% 
  arrange(chromosome, 
          position)

gene_name_1 = SKR_gene_overlap_tib %>% 
  # pull(gene_id) %>% 
  as_tibble() %>% 
  dplyr::select(chromosome, 
                position, 
                start,
                i.start, 
                end, 
                i.end,
                gene_id,
                feature,
                FST) %>% 
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
                i.start,
                end, 
                i.end,
                FST,
                feature,
                gene_name) %>% 
  na.omit()

gene_name_2 = SKR_gene_overlap_tib %>% 
  # pull(gene_id) %>% 
  as_tibble() %>% 
  dplyr::select(chromosome, 
                position, 
                start,
                i.start,
                end, 
                i.end,
                gene_id,
                feature,
                FST,) %>% 
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
                i.start,
                end, 
                i.end,
                FST,
                feature,
                gene_name) %>% 
  na.omit()


SKR_FST_out_genes = bind_rows(gene_name_1, 
                               gene_name_2) %>% 
  arrange(chromosome, 
          position) %>% 
  distinct(gene_name, 
           .keep_all = T) %>% 
  filter(!grepl('ENSG', 
                gene_name)) %>% 
  filter(feature == 'gene')

SKR_FST_out_genes %>% 
  write_csv('SKR_FST_0.5%_outlier_genes.csv')

SKR_FST_out_genes %>% 
  select(gene_name) %>% 
  write_tsv('SKR_FST_0.5%_outlier_gene_names_only.tsv', 
            col_names = F)


#### GTS_CSWY outliers ----------------------------------------------------

GTS_CSWY_div_snps = read_csv('GTS_CSWY_TOP_DAWG_Fst_clean.csv') %>% 
  select(-NMISS, 
         -FST) %>%
  rename(position = POS, 
         chromosome = CHR, 
         FST = FST_zero) %>% 
  filter(value == 'Outlier') %>% 
  arrange(chromosome, 
          position) %>% 
  group_by(chromosome)
## create a 100bp window around the snp identified
## settting the start and end range for the methylation data
## this number can be as big or small as you want depending
## on how liberal you want it

GTS_CSWY_div_window = GTS_CSWY_div_snps %>%
  group_by(chromosome)%>%
  mutate(start = position-100,
         end = position+100) %>% 
  separate(col = chromosome, 
           into = c('chr', 
                    'chr_name'), 
           sep = '_') %>% 
  unite(chromosome, 
        chr:chr_name,
        sep = '',
        remove = F) %>% 
  select(-chr, 
         -chr_name, 
         -SNP) %>% 
  group_by(chromosome)

setDT(GTS_CSWY_div_window)
setDT(gene_annotation)

setkey(GTS_CSWY_div_window, 
       chromosome, 
       start, 
       end)

GTS_CSWY_gene_overlap = foverlaps(gene_annotation,
                              GTS_CSWY_div_window,
                              # by.x = start,
                              # by.y = end,
                              type="any")


GTS_CSWY_gene_overlap_tib = as_tibble(GTS_CSWY_gene_overlap) %>% 
  na.omit() %>% 
  filter(chromosome != 'chrUn') %>% 
  arrange(chromosome, 
          position)

gene_name_1 = GTS_CSWY_gene_overlap_tib %>% 
  # pull(gene_id) %>% 
  as_tibble() %>% 
  dplyr::select(chromosome, 
                position, 
                start,
                i.start, 
                end, 
                i.end,
                gene_id,
                feature,
                FST) %>% 
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
                i.start,
                end, 
                i.end,
                FST,
                feature,
                gene_name) %>% 
  na.omit()

gene_name_2 = GTS_CSWY_gene_overlap_tib %>% 
  # pull(gene_id) %>% 
  as_tibble() %>% 
  dplyr::select(chromosome, 
                position, 
                start,
                i.start,
                end, 
                i.end,
                gene_id,
                feature,
                FST,) %>% 
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
                i.start,
                end, 
                i.end,
                FST,
                feature,
                gene_name) %>% 
  na.omit()


GTS_CSWY_FST_out_genes = bind_rows(gene_name_1, 
                               gene_name_2) %>% 
  arrange(chromosome, 
          position) %>% 
  distinct(gene_name, 
           .keep_all = T) %>% 
  filter(!grepl('ENSG', 
                gene_name)) %>% 
  filter(feature == 'gene')

GTS_CSWY_FST_out_genes %>% 
  write_csv('GTS_CSWY_FST_0.5%_outlier_genes.csv')


GTS_CSWY_FST_out_genes %>% 
  select(gene_name) %>% 
  write_tsv('GTS_CSWY_FST_0.5%_outlier_gene_names_only.tsv', 
            col_names = F)






# Gene overlap b/w pops ---------------------------------------------------

ASHN_fst_genes = read_csv('ASHN_FST_0.5%_outlier_genes.csv') %>%
  select(gene_name)
MYV_fst_genes = read_csv('MYV_FST_0.5%_outlier_genes.csv') %>% 
  select(gene_name)
SKR_fst_genes = read_csv('SKR_FST_0.5%_outlier_genes.csv') %>% 
  select(gene_name)
GTS_CSWY_fst_genes =  read_csv('GTS_CSWY_FST_0.5%_outlier_genes.csv') %>% 
  select(gene_name)


intersect(ASHN_fst_genes, 
          MYV_fst_genes)

intersect(ASHN_fst_genes, 
          SKR_fst_genes)
## a single overlapping gene

intersect(ASHN_fst_genes, 
          GTS_CSWY_fst_genes)
##single gene overlapping

intersect(MYV_fst_genes, 
          SKR_fst_genes)
## 4 overlapping genes

intersect(MYV_fst_genes, 
          GTS_CSWY_fst_genes)
## 3 genes overlapping

intersect(SKR_fst_genes, 
          GTS_CSWY_fst_genes)
## 3 overlapping genes


# Gene functions b/w pops -------------------------------------------------

ASHN_fst_genes = read_csv('ASHN_FST_0.5%_outlier_genes.csv') 
MYV_fst_genes = read_csv('MYV_FST_0.5%_outlier_genes.csv') 
SKR_fst_genes = read_csv('SKR_FST_0.5%_outlier_genes.csv') 
GTS_CSWY_fst_genes = read_csv('GTS_CSWY_FST_0.5%_outlier_genes.csv') 


ASHN_fst_genes %>% 
  select(chromosome, 
         position, 
         feature, 
         gene_name) %>% 
  group_by(feature) %>% 
  summarize(n_feature = n())

MYV_fst_genes %>% 
  select(chromosome, 
         position, 
         feature, 
         gene_name) %>% 
  group_by(feature) %>% 
  summarize(n_feature = n())

SKR_fst_genes %>% 
  select(chromosome, 
         position, 
         feature, 
         gene_name) %>% 
  group_by(feature) %>% 
  summarize(n_feature = n())

GTS_CSWY_fst_genes %>% 
  select(chromosome, 
         position, 
         feature, 
         gene_name) %>% 
  group_by(feature) %>% 
  summarize(n_feature = n())


# ASHN genes no 100bp window ------------------------------

ASHN_div_window = ASHN_div_snps %>%
  group_by(chromosome)%>%
  mutate(start = position-1,
         end = position+1) %>% 
  separate(col = chromosome, 
           into = c('chr', 
                    'chr_name'), 
           sep = '_') %>% 
  unite(chromosome, 
        chr:chr_name,
        sep = '',
        remove = F) %>% 
  select(-chr, 
         -chr_name, 
         -SNP) %>% 
  group_by(chromosome)


setDT(ASHN_div_window)
setDT(gene_annotation)

setkey(ASHN_div_window, 
       chromosome, 
       start, 
       end)

ASHN_gene_overlap = foverlaps(gene_annotation,
                              ASHN_div_window,
                              # by.x = start,
                              # by.y = end,
                              type="any")


ASHN_gene_overlap_tib = as_tibble(ASHN_gene_overlap) %>% 
  na.omit() %>% 
  filter(chromosome != 'chrUn') %>% 
  arrange(chromosome, 
          position)

ASHN_gene_overlap_tib %>% 
  filter(feature == 'CDS') %>% 
  pull(gene_id)

ASHN_gene_overlap_tib %>% 
  group_by(feature) %>% 
  summarize(n = n())

gene_name_1 = ASHN_gene_overlap_tib %>% 
  # pull(gene_id) %>% 
  as_tibble() %>% 
  dplyr::select(chromosome, 
                position, 
                start,
                i.start, 
                end, 
                i.end,
                gene_id,
                feature,
                FST) %>% 
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
                i.start,
                end, 
                i.end,
                FST,
                feature,
                gene_name) %>% 
  na.omit()

ASHN_gene_table = gene_name_1 %>%
  arrange(chromosome, 
          position) %>% 
  distinct(gene_name, 
           .keep_all = T) %>%
  group_by(feature) %>%
  summarize(n = n())

# gene_name_1 %>%
#   arrange(chromosome, 
#           position) %>% 
#   distinct(gene_name, 
#            .keep_all = T) %>%
#   filter(feature %in% c('gene',
#                         'CDS')) %>%
#   arrange(chromosome,
#           position) %>%
#   # filter(feature == 'gene') %>%
#   # filter(!grepl('ENSG',
#   #               gene_name))
#   select(gene_name)%>% 
#   write_csv('ASHN_NoWindow_FST_0.5%_outlier_genes.csv')


ASHN_gene_only = gene_name_1 %>%
  arrange(chromosome, 
          position) %>% 
  distinct(gene_name, 
           .keep_all = T) %>% 
  filter(feature == 'gene') %>% 
    # filter(!grepl('ENSG',
    #               gene_name))
  select(gene_name)%>% 
  write_csv('ASHN_NoWindow_FST_0.5%_outlier_genes.csv')

ASHN_regulatory_coding_genes = gene_name_1 %>%
  arrange(chromosome, 
          position) %>% 
  distinct(gene_name, 
           .keep_all = T) %>% 
  select(gene_name)%>% 
  write_csv('ASHN_NoWindow_FST_0.5%_outlier_regulatory_and_genes.csv')

# 
#   filter(!grepl('ENSG', 
#                 gene_name)) 


# gene_name_2 = ASHN_gene_overlap_tib %>% 
#   # pull(gene_id) %>% 
#   as_tibble() %>% 
#   dplyr::select(chromosome, 
#                 position, 
#                 start,
#                 i.start,
#                 end, 
#                 i.end,
#                 gene_id,
#                 feature,
#                 FST,) %>% 
#   separate(col = gene_id, 
#            into = c('ensemble_id', 
#                     'gene_name', 
#                     'parent_code', 
#                     'gene_name2'), 
#            sep = ';') %>% 
#   separate(col = parent_code, 
#            into = c('Garbage', 
#                     'gene_name'), 
#            sep = '=') %>% 
#   dplyr::select(chromosome, 
#                 position, 
#                 start, 
#                 i.start,
#                 end, 
#                 i.end,
#                 FST,
#                 feature,
#                 gene_name) %>% 
#   na.omit()
# 
# gene_name_2 %>% 
#   group_by(feature) %>% 
#   summarize(n = n())

# ASHN_FST_out_genes = bind_rows(gene_name_1, 
#                                gene_name_2) %>% 
#   arrange(chromosome, 
#           position) %>% 
#   distinct(gene_name, 
#            .keep_all = T) %>% 
#   filter(!grepl('ENSG', 
#                 gene_name)) 
# 
# ASHN_FST_out_genes %>% 
#   group_by(feature) %>% 
#   summarize(n_features = n())

# ASHN_FST_out_genes %>% 
#     filter(feature == 'gene') %>% 
#   write_csv('ASHN_NoWindow_FST_0.5%_outlier_genes.csv')
# 
# ASHN_FST_out_genes %>% 
#   select(gene_name) %>% 
#   write_tsv('ASHN_NoWindow_FST_0.5%_outlier_gene_names_only.tsv', 
#             col_names = F)



##
# MYV genes no 100bp window ------------------------------


MYV_div_snps = read_csv('MYV_TOP_DAWG_Fst_clean.csv') %>% 
  select(-NMISS, 
         -FST) %>%
  rename(position = POS, 
         chromosome = CHR, 
         FST = FST_zero) %>% 
  filter(value == 'Outlier') %>% 
  arrange(chromosome, 
          position) %>% 
  group_by(chromosome)

MYV_div_window = MYV_div_snps %>%
  group_by(chromosome)%>%
  mutate(start = position-1,
         end = position+1) %>% 
  separate(col = chromosome, 
           into = c('chr', 
                    'chr_name'), 
           sep = '_') %>% 
  unite(chromosome, 
        chr:chr_name,
        sep = '',
        remove = F) %>% 
  select(-chr, 
         -chr_name, 
         -SNP) %>% 
  group_by(chromosome)


setDT(MYV_div_window)
setDT(gene_annotation)

setkey(MYV_div_window, 
       chromosome, 
       start, 
       end)

MYV_gene_overlap = foverlaps(gene_annotation,
                              MYV_div_window,
                              # by.x = start,
                              # by.y = end,
                              type="any")


MYV_gene_overlap_tib = as_tibble(MYV_gene_overlap) %>% 
  na.omit() %>% 
  filter(chromosome != 'chrUn') %>% 
  arrange(chromosome, 
          position)

gene_name_1 = MYV_gene_overlap_tib %>% 
  # pull(gene_id) %>% 
  as_tibble() %>% 
  dplyr::select(chromosome, 
                position, 
                start,
                i.start, 
                end, 
                i.end,
                gene_id,
                feature,
                FST) %>% 
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
                i.start,
                end, 
                i.end,
                FST,
                feature,
                gene_name) %>% 
  na.omit()

MYV_gene_table = gene_name_1 %>%
  arrange(chromosome, 
          position) %>% 
  distinct(gene_name, 
           .keep_all = T) %>%
  group_by(feature) %>%
  summarize(n = n())

gene_name_1 %>%
  arrange(chromosome, 
          position) %>% 
  distinct(gene_name, 
           .keep_all = T) %>%
  filter(feature %in% c('gene',
                        'CDS')) %>%
  arrange(chromosome,
          position) %>% View()
  # filter(feature == 'gene') %>%
  # filter(!grepl('ENSG',
  #               gene_name))
  select(gene_name)%>% 
  write_csv('ASHN_NoWindow_FST_0.5%_outlier_genes.csv')


MYV_gene_only = gene_name_1 %>%
  arrange(chromosome, 
          position) %>% 
  distinct(gene_name, 
           .keep_all = T) %>% 
  filter(feature == 'gene') %>% 
  select(gene_name)%>% 
  write_csv('MYV_NoWindow_FST_0.5%_outlier_genes.csv')

MYV_regulatory_coding_genes = gene_name_1 %>%
  arrange(chromosome, 
          position) %>% 
  distinct(gene_name, 
           .keep_all = T) %>% 
  select(gene_name)%>% 
  write_csv('MYV_NoWindow_FST_0.5%_outlier_regulatory_and_genes.csv')

# 
# 
# 
# gene_name_2 = MYV_gene_overlap_tib %>% 
#   # pull(gene_id) %>% 
#   as_tibble() %>% 
#   dplyr::select(chromosome, 
#                 position, 
#                 start,
#                 i.start,
#                 end, 
#                 i.end,
#                 gene_id,
#                 feature,
#                 FST,) %>% 
#   separate(col = gene_id, 
#            into = c('ensemble_id', 
#                     'gene_name', 
#                     'parent_code', 
#                     'gene_name2'), 
#            sep = ';') %>%
#   separate(col = parent_code, 
#            into = c('Garbage', 
#                     'gene_name'), 
#            sep = '=') %>% 
#   dplyr::select(chromosome, 
#                 position, 
#                 start, 
#                 i.start,
#                 end, 
#                 i.end,
#                 FST,
#                 feature,
#                 gene_name) %>% 
#   na.omit()
# 
# 
# MYV_FST_out_genes = bind_rows(gene_name_1, 
#                                gene_name_2) %>% 
#   arrange(chromosome, 
#           position) %>% 
#   distinct(gene_name, 
#            .keep_all = T) %>% 
#   filter(!grepl('ENSG', 
#                 gene_name)) 
# 
# MYV_FST_out_genes %>% 
#   group_by(feature) %>% 
#   summarize(n_features = n())
# 
# MYV_FST_out_genes %>% 
#   filter(feature == 'gene') %>% 
#   write_csv('MYV_NoWindow_FST_0.5%_outlier_genes.csv')
# 
# MYV_FST_out_genes %>% 
#   select(gene_name) %>% 
#   write_tsv('MYV_NoWindow_FST_0.5%_outlier_gene_names_only.tsv', 
#             col_names = F)

# SKR genes no 100bp window ------------------------------
SKR_div_snps = read_csv('SKR_TOP_DAWG_Fst_clean.csv') %>% 
  select(-NMISS, 
         -FST) %>%
  rename(position = POS, 
         chromosome = CHR, 
         FST = FST_zero) %>% 
  filter(value == 'Outlier') %>% 
  arrange(chromosome, 
          position) %>% 
  group_by(chromosome)

SKR_div_window = SKR_div_snps %>%
  group_by(chromosome)%>%
  mutate(start = position-1,
         end = position+1) %>% 
  separate(col = chromosome, 
           into = c('chr', 
                    'chr_name'), 
           sep = '_') %>% 
  unite(chromosome, 
        chr:chr_name,
        sep = '',
        remove = F) %>% 
  select(-chr, 
         -chr_name, 
         -SNP) %>% 
  group_by(chromosome)


setDT(SKR_div_window)
setDT(gene_annotation)

setkey(SKR_div_window, 
       chromosome, 
       start, 
       end)

SKR_gene_overlap = foverlaps(gene_annotation,
                              SKR_div_window,
                              # by.x = start,
                              # by.y = end,
                              type="any")


SKR_gene_overlap_tib = as_tibble(SKR_gene_overlap) %>% 
  na.omit() %>% 
  filter(chromosome != 'chrUn') %>% 
  arrange(chromosome, 
          position)

gene_name_1 = SKR_gene_overlap_tib %>% 
  # pull(gene_id) %>% 
  as_tibble() %>% 
  dplyr::select(chromosome, 
                position, 
                start,
                i.start, 
                end, 
                i.end,
                gene_id,
                feature,
                FST) %>% 
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
                i.start,
                end, 
                i.end,
                FST,
                feature,
                gene_name) %>% 
  na.omit()




SKR_gene_table = gene_name_1 %>%
  arrange(chromosome, 
          position) %>% 
  distinct(gene_name, 
           .keep_all = T) %>%
  group_by(feature) %>%
  summarize(n = n())

gene_name_1 %>%
  arrange(chromosome, 
          position) %>% 
  distinct(gene_name, 
           .keep_all = T) %>%
  filter(feature %in% c('gene',
                        'CDS')) %>%
  arrange(chromosome,
          position) %>% View()

SKR_gene_only = gene_name_1 %>%
  arrange(chromosome, 
          position) %>% 
  distinct(gene_name, 
           .keep_all = T) %>% 
  filter(feature == 'gene') %>% 
  select(gene_name)%>% 
  write_csv('SKR_NoWindow_FST_0.5%_outlier_genes.csv')

SKR_regulatory_coding_genes = gene_name_1 %>%
  arrange(chromosome, 
          position) %>% 
  distinct(gene_name, 
           .keep_all = T) %>% 
  select(gene_name)%>% 
  write_csv('SKR_NoWindow_FST_0.5%_outlier_regulatory_and_genes.csv')

# 

# gene_name_2 = SKR_gene_overlap_tib %>% 
#   # pull(gene_id) %>% 
#   as_tibble() %>% 
#   dplyr::select(chromosome, 
#                 position, 
#                 start,
#                 i.start,
#                 end, 
#                 i.end,
#                 gene_id,
#                 feature,
#                 FST,) %>% 
#   separate(col = gene_id, 
#            into = c('ensemble_id', 
#                     'gene_name', 
#                     'parent_code', 
#                     'gene_name2'), 
#            sep = ';') %>%
#   separate(col = parent_code, 
#            into = c('Garbage', 
#                     'gene_name'), 
#            sep = '=') %>% 
#   dplyr::select(chromosome, 
#                 position, 
#                 start, 
#                 i.start,
#                 end, 
#                 i.end,
#                 FST,
#                 feature,
#                 gene_name) %>% 
#   na.omit()
# 
# 
# SKR_FST_out_genes = bind_rows(gene_name_1, 
#                                gene_name_2) %>% 
#   arrange(chromosome, 
#           position) %>% 
#   distinct(gene_name, 
#            .keep_all = T) %>% 
#   filter(!grepl('ENSG', 
#                 gene_name)) 
# 
# SKR_FST_out_genes %>% 
#   group_by(feature) %>% 
#   summarize(n_features = n())
# 
# SKR_FST_out_genes %>% 
#   filter(feature == 'gene') %>% 
#   write_csv('SKR_NoWindow_FST_0.5%_outlier_genes.csv')
# 
# SKR_FST_out_genes %>% 
#   select(gene_name) %>% 
#   write_tsv('SKR_NoWindow_FST_0.5%_outlier_gene_names_only.tsv', 
#             col_names = F)

# GTS_CSWY genes no 100bp window ------------------------------
GTS_CSWY_div_snps = read_csv('GTS_CSWY_TOP_DAWG_Fst_clean.csv') %>% 
  select(-NMISS, 
         -FST) %>%
  rename(position = POS, 
         chromosome = CHR, 
         FST = FST_zero) %>% 
  filter(value == 'Outlier') %>% 
  arrange(chromosome, 
          position) %>% 
  group_by(chromosome)

GTS_CSWY_div_window = GTS_CSWY_div_snps %>%
  group_by(chromosome)%>%
  mutate(start = position-1,
         end = position+1) %>% 
  separate(col = chromosome, 
           into = c('chr', 
                    'chr_name'), 
           sep = '_') %>% 
  unite(chromosome, 
        chr:chr_name,
        sep = '',
        remove = F) %>% 
  select(-chr, 
         -chr_name, 
         -SNP) %>% 
  group_by(chromosome)


setDT(GTS_CSWY_div_window)
setDT(gene_annotation)

setkey(GTS_CSWY_div_window, 
       chromosome, 
       start, 
       end)

GTS_CSWY_gene_overlap = foverlaps(gene_annotation,
                              GTS_CSWY_div_window,
                              # by.x = start,
                              # by.y = end,
                              type="any")


GTS_CSWY_gene_overlap_tib = as_tibble(GTS_CSWY_gene_overlap) %>% 
  na.omit() %>% 
  filter(chromosome != 'chrUn') %>% 
  arrange(chromosome, 
          position)

gene_name_1 = GTS_CSWY_gene_overlap_tib %>% 
  # pull(gene_id) %>% 
  as_tibble() %>% 
  dplyr::select(chromosome, 
                position, 
                start,
                i.start, 
                end, 
                i.end,
                gene_id,
                feature,
                FST) %>% 
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
                i.start,
                end, 
                i.end,
                FST,
                feature,
                gene_name) %>% 
  na.omit()

GTS_GAR_gene_table = gene_name_1 %>%
  arrange(chromosome, 
          position) %>% 
  distinct(gene_name, 
           .keep_all = T) %>%
  group_by(feature) %>%
  summarize(n = n())

gene_name_1 %>%
  arrange(chromosome, 
          position) %>% 
  distinct(gene_name, 
           .keep_all = T) %>%
  filter(feature %in% c('gene',
                        'CDS')) %>%
  arrange(chromosome,
          position) %>% View()

GTS_GAR_gene_only = gene_name_1 %>%
  arrange(chromosome, 
          position) %>% 
  distinct(gene_name, 
           .keep_all = T) %>% 
  filter(feature == 'gene') %>% 
  select(gene_name)%>% 
  write_csv('GTS_GAR_NoWindow_FST_0.5%_outlier_genes.csv')

GTS_GAR_regulatory_coding_genes = gene_name_1 %>%
  arrange(chromosome, 
          position) %>% 
  distinct(gene_name, 
           .keep_all = T) %>% 
  select(gene_name)%>% 
  write_csv('GTS_GAR_NoWindow_FST_0.5%_outlier_regulatory_and_genes.csv')

# 
# gene_name_2 = GTS_CSWY_gene_overlap_tib %>% 
#   # pull(gene_id) %>% 
#   as_tibble() %>% 
#   dplyr::select(chromosome, 
#                 position, 
#                 start,
#                 i.start,
#                 end, 
#                 i.end,
#                 gene_id,
#                 feature,
#                 FST,) %>% 
#   separate(col = gene_id, 
#            into = c('ensemble_id', 
#                     'gene_name', 
#                     'parent_code', 
#                     'gene_name2'), 
#            sep = ';') %>%
#   separate(col = parent_code, 
#            into = c('Garbage', 
#                     'gene_name'), 
#            sep = '=') %>% 
#   dplyr::select(chromosome, 
#                 position, 
#                 start, 
#                 i.start,
#                 end, 
#                 i.end,
#                 FST,
#                 feature,
#                 gene_name) %>% 
#   na.omit()
# 
# 
# GTS_CSWY_FST_out_genes = bind_rows(gene_name_1, 
#                                gene_name_2) %>% 
#   arrange(chromosome, 
#           position) %>% 
#   distinct(gene_name, 
#            .keep_all = T) %>% 
#   filter(!grepl('ENSG', 
#                 gene_name)) 
# 
# GTS_CSWY_FST_out_genes %>% 
#   group_by(feature) %>% 
#   summarize(n_features = n())
# 
# GTS_CSWY_FST_out_genes %>% 
#   filter(feature == 'gene') %>% 
#   write_csv('GTS_CSWY_NoWindow_FST_0.5%_outlier_genes.csv')
# 
# GTS_CSWY_FST_out_genes %>% 
#   select(gene_name) %>% 
#   write_tsv('GTS_CSWY_NoWindow_FST_0.5%_outlier_gene_names_only.tsv', 
#             col_names = F)



# Gene overlap b/w pops ---------------------------------------------------

ASHN_genes = read_csv('ASHN_NoWindow_FST_0.5%_outlier_genes.csv') %>% 
  select(gene_name)
MYV_genes = read_csv('MYV_NoWindow_FST_0.5%_outlier_genes.csv') %>% 
  select(gene_name)
SKR_genes = read_csv('SKR_NoWindow_FST_0.5%_outlier_genes.csv') %>% 
  select(gene_name)
GTS_CSWY_genes = read_csv('GTS_Gar_NoWindow_FST_0.5%_outlier_genes.csv') %>% 
  select(gene_name)

intersect(ASHN_genes, 
          MYV_genes)

intersect(ASHN_genes, 
          SKR_genes)

intersect(ASHN_genes, 
          GTS_CSWY_genes)

intersect(MYV_genes, 
          SKR_genes)

intersect(MYV_genes, 
          GTS_CSWY_genes)

intersect(SKR_genes, 
          GTS_CSWY_genes)


##
# Methylation outliers ----------------------------------------------------


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

View(gene_overlap_tib)

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
