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
         skip = 1) %>% 
  select(X3) %>% 
  distinct()


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

genes = gene_annotation %>% 
  arrange(chromosome) %>% 
  filter(feature == 'gene') 
# %>% 
#   select(chromosome, 
#          feature, 
#          gene_id, 
#          position)

ASHN_div_snps = read_csv('ASHN_TOP_DAWG_Fst_clean.csv') %>% 
  select(-NMISS, 
         -FST) %>%
  rename(position = POS, 
         chromosome = CHR, 
         FST = FST_zero) %>% 
  filter(value == 'Outlier') %>% 
  group_by(chromosome)

ASHN_div_snps %>% 
  summarise(snp_per_chr = n()) %>% 
  View()

## Cannot figure this out right now
## Need to find out if the average positions of each gene
## lies within the window range
## my brain hurts

inner_join(WC_Fst_regions, 
           genes, 
           by = c('chromosome', 
                  'start', 
                  'end'))


gene_test = genes %>% 
  filter(chromosome == 'chrI')


test_data = WC_Fst_regions %>% 
  filter(chromosome == 'chr_I')


# Test old R functions  ---------------------------------------------------

library(tidyverse)
#its the only package from the original list we actually need >:(

#alter this line to the directory on your computer containing the files to be analyzed

# this function takes a .gff or .gff3 file and creates a tibble
# arguments:
# filename -  the name of the gff file
# lead_skip - this argument specifies the
# chromosome - the chromosome to subset from the full gff, default is 'all' and no subset is performed.
read_gff = function(filename, lead_skip = 9, chromosome = "all"){
  #CAM - personally I would pass a vector of names to the read_tsv as opposed to renaming in the pipe (this could be an argument too).
  colnames = c("chr", "source", "feature", "start", "end", "score", "strand", "frame", "attribute" )
  #read in the raw file with readr
  gff = read_tsv(filename, col_names = colnames, skip = lead_skip) %>%
    filter(feature == 'gene') %>% 
    arrange(start, end) %>% 
    mutate(mid = start + (end-start)/2) ## arrange each gene by its start and end points on each chromosome
  #subset if user changed the chromosome argument
  if(chromosome != "all"){
    sub_gff = gff %>% filter(chr == chromosome)
    return(sub_gff)
  }
  #if the default chr arg was not changed, return full gff
  return(gff)
}

# this function takes a tab delimited file of outlier loci and creates a tibble
# arguments:
# filename -  the name of the gff file
# chromosome - the chromosome to subset from the full gff, default is 'all' and no subset is performed.
read_outlier = function(filename, chromosome = "all"){
  #read in the dataframe
  out_dat = read_csv(filename) %>% 
    arrange(chromosome) 
  # %>% 
  #   filter(chromo != 'Contigs')
  # #subset if user changed the chromosome argument
  if(chromosome != "all"){
    sub_out = out_dat %>% 
      filter(chromosome == chromosome) 
    return(sub_out)
  }
  #if the default chr arg was not changed, return full gff
  return(out_dat)
}



#I set important filenames etc. to variables, so they can be easily found, changed, and reused.
gff_filename = "stickleback_v5_ensembl_genes.gff3.gz"
ncbi_code_AC08 = "chr_I"
outlier_data = "ASHN_TOP_DAWG_FST_outlier.csv"
outlier_code_AC08 = "chr_I"
# read_csv('ASHN_TOP_DAWG_Fst_clean.csv') %>% 
#   select(-NMISS, 
#          -FST) %>%
#   rename(position = POS, 
#          chromosome = CHR, 
#          FST = FST_zero) %>% 
#   filter(value == 'Outlier') %>% 
#   write_csv('ASHN_TOP_DAWG_FST_outlier.csv')

########
# read in the data
# all the detail above is abstracted away into these two lines
gff = read_gff(gff_filename, chromosome = 'all') %>% 
  filter(chr != 'chrM') %>% 
  arrange(chr, 
          mid) %>% 
  group_by(chr)
outliers = read_outlier(outlier_data, chromosome = 'all') %>% 
  arrange(chromosome, 
          position) %>% 
  separate(col = chromosome, 
           into = c('chr', 'chr_name'), 
           sep = '_') %>% 
  unite(chromosome, 
        chr:chr_name,
        sep = '',
        remove = F) %>% 
  select(-chr, 
         -chr_name) %>% 
  group_by(chromosome)

## might have to use a map function to see if these actually line up



# obtain outlier positions
pos = outliers$position



#obtain genes in those positions
gene_regions = gff %>% 
  group_by(chr) %>% 
  mutate(hit_dist = abs(mid - pos)) %>% 
  arrange(hit_dist) %>% 
  filter(hit_dist < 1000) %>%
  pull(attribute)
  # select(chr, start, end, attribute, hit_dist) %>% 
  # pull(attribute)


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
