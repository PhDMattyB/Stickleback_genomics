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

## Reads in the gff file that has the stickleback gene annotations
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

## Take a look at what the functional annotations actually are
## CDS is coding sequence, this is the actual coding region
## gene means it's in the gene body
## mRNA means it's associated with an mRNA transcript
## three_prime_UTR means it's in the non-coding tail region
## five_primte_UTR is the promoter

gene_annotation %>% 
  arrange(chromosome) %>%
  # filter(feature == 'gene')
  group_by(feature) %>% 
  distinct(feature)

## Next read in the outlier loci you have
## can find genes within a window or just at the site of interest. 

## Site of interest

ASHN_div_snps = read_csv('ASHN_TOP_DAWG_Fst_clean.csv') %>% 
  select(-NMISS, 
         -FST) %>%
  rename(position = POS, 
         chromosome = CHR, 
         FST = FST_zero) %>% 
  filter(value == 'Outlier') %>% 
  arrange(chromosome, 
          position) %>% 
  group_by(chromosome)%>%
  mutate(start = position-1,   ## can change this to whatever window of interest you want around the site of interest
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


## using data.table to set key value pair
setDT(ASHN_div_window)
setDT(gene_annotation)

setkey(ASHN_div_window, 
       chromosome, 
       start, 
       end)

## aligns the sites of interest with the annotated genome
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

  gene_name = ASHN_gene_overlap_tib %>% 
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

ASHN_gene_table = gene_name %>%
  arrange(chromosome, 
          position) %>% 
  distinct(gene_name, 
           .keep_all = T) %>%
  group_by(feature) %>%
  summarize(n = n())


## can modify the script below to pull out any feature you want
  ASHN_gene_only = gene_name %>%
  arrange(chromosome, 
          position) %>% 
  distinct(gene_name, 
           .keep_all = T) %>% 
  filter(feature == 'gene') %>% 
    filter(!grepl('ENSG',
                  gene_name)) %>% 
  select(gene_name)%>% 
  write_csv('ASHN_NoWindow_FST_0.5%_outlier_genes.csv')

ASHN_regulatory_coding_genes = gene_name %>%
  arrange(chromosome, 
          position) %>% 
  distinct(gene_name, 
           .keep_all = T) %>% 
  select(gene_name)%>% 
  write_csv('ASHN_NoWindow_FST_0.5%_outlier_regulatory_and_genes.csv')


