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
library(enrichR)
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

all_gene_names = gene_annotation %>% 
  # filter(chromosome == 'chrXXI') %>% 
  separate(col = gene_id, 
           into = c('ensemble_id', 
                    'gene_name', 
                    'parent_code'), 
           sep = ';') %>% 
  separate(col = ensemble_id, 
           into = c('garbage', 
                    'ensemble_id'), 
           sep = '=') %>% 
  dplyr::select(-garbage) %>% 
  separate(col = gene_name, 
           into = c('Garbage', 
                    'gene_name'), 
           sep = '=') %>% 
  dplyr::select(-Garbage) %>% 
  separate(col = parent_code, 
           into = c('garbage', 
                    'parent_gene_name'), 
           sep = '=') %>% 
  dplyr::select(-garbage) %>% 
  filter(feature == 'gene') %>% 
  ungroup(chromosome) %>% 
  dplyr::select(gene_name) %>% 
  distinct(gene_name) %>% 
  filter(!if_any(everything(), ~ grepl('ENSG', .))) %>% 
  mutate_all(as.character) %>% 
  as.data.frame()


# ChrXXI inversion genes --------------------------------------------------

gene_annotation %>% 
  filter(chromosome == 'chrXXI') %>% 
# %>%
#   separate(col = gene_id, 
#            into = c('ensemble_id', 
#                     'gene_name', 
#                     'parent_code'), 
#            sep = ';') %>% 
#   separate(col = ensemble_id, 
#            into = c('garbage', 
#                     'ensemble_id'), 
#            sep = '=') %>% 
#   dplyr::select(-garbage) %>% 
#   separate(col = gene_name, 
#            into = c('Garbage', 
#                     'gene_name'), 
#            sep = '=') %>% 
#   dplyr::select(-Garbage) %>% 
#   separate(col = parent_code, 
#            into = c('garbage', 
#                     'parent_gene_name'), 
#            sep = '=') %>% 
#   dplyr::select(-garbage) %>% 
  filter(start >= 9963830) %>% 
  filter(end <= 11574024)

chrxxi_inversion_genes = gene_annotation %>% 
  filter(chromosome == 'chrXXI') %>% 
  filter(start >= 9963830) %>% 
  filter(end <= 11574024) %>% 
  separate(col = gene_id, 
           into = c('ensemble_id', 
                    'gene_name', 
                    'parent_code'), 
           sep = ';') %>% 
  separate(col = ensemble_id, 
           into = c('garbage', 
                    'ensemble_id'), 
           sep = '=') %>% 
  dplyr::select(-garbage) %>% 
  separate(col = gene_name, 
           into = c('Garbage', 
                    'gene_name'), 
           sep = '=') %>% 
  dplyr::select(-Garbage) %>% 
  separate(col = parent_code, 
           into = c('garbage', 
                    'parent_gene_name'), 
           sep = '=') %>% 
  dplyr::select(-garbage) 

chrxxi_inversion_genes %>% 
  # filter(feature == 'CDS') %>% View()
  # filter(feature %in% c('gene', 'CDS')) %>% View()
  filter(feature == 'gene') %>%
  ungroup() %>% 
  select(gene_name) %>% 
  write_tsv('CHRXXI_Inversion_gene_names_fixed.txt')
  
distinct(gene_name) %>% 
  filter(!if_any(everything(), ~ grepl('ENSG', .)))


# chr xxi inversion gene ontology -----------------------------------------
# library(enrichR)
listEnrichrSites()
setEnrichrSite("Enrichr")

database_list <- listEnrichrDbs()

database_list %>% 
  as_tibble() %>% 
  View()
#GO_Biological_Process_2025
#GO_Cellular_Component_2025
#	GO_Molecular_Function_2025
#	KEGG_2021_Human
#PhenGenI_Association_2021
#Reactome_Pathways_2024
#TRANSFAC_and_JASPAR_PWMs
#WikiPathways_2024_Mouse
#	WikiPathways_2024_Human

dbs <- c("GO_Molecular_Function_2025", 
         "GO_Cellular_Component_2025",
         "GO_Biological_Process_2025",
         'KEGG_2021_Human', 
         'PhenGenI_Association_2021', 
         'Reactome_Pathways_2024', 
         'TRANSFAC_and_JASPAR_PWMs', 
         'WikiPathways_2024_Mouse', 
         'WikiPathways_2024_Human')

dbs2 = c('PerturbAtlas')

chrxxi_inversion_gene_names = read_tsv('CHRXXI_Inversion_gene_names_fixed.txt')%>%
  # as_tibble() %>% 
  mutate_all(as.character) %>% 
  as.data.frame()

# data(input)
enriched = enrichR::enrichr(chrxxi_inversion_gene_names, 
                    dbs)
head(enriched[["GO_Biological_Process_2025"]])
enriched$GO_Cellular_Component_2025 %>% 
  as_tibble() %>% 
  filter(Adjusted.P.value <= 0.05)
enriched$GO_Molecular_Function_2025 %>% 
  as_tibble() %>% 
  filter(Adjusted.P.value <= 0.05)
enriched$GO_Biological_Process_2025 %>% 
  as_tibble() %>% 
  filter(Adjusted.P.value <= 0.05)
enriched$KEGG_2021_Human %>% 
  as_tibble() %>% 
  filter(Adjusted.P.value <= 0.05)
enriched$PhenGenI_Association_2021 %>% 
  as_tibble() %>% 
  filter(Adjusted.P.value <= 0.05)
enriched$Reactome_Pathways_2024 %>% 
  as_tibble() %>% 
  filter(Adjusted.P.value <= 0.05)
enriched$WikiPathways_2024_Human %>% 
  as_tibble() %>% 
  filter(Adjusted.P.value <= 0.05)

chrxxi_background = pull(chrxxi_inversion_gene_names, 
                         gene_name)
stickle_background = pull(all_gene_names, 
                         gene_name)

enriched_chrxxi_inversion <- enrichr(chrxxi_inversion_gene_names, 
                     dbs, 
                     background = stickle_background, 
                     include_overlap = T)

# GO_mol_func = enriched_chrxxi_inversion$GO_Molecular_Function_2025%>% 
#   as_tibble() 
# p.adjust(GO_mol_func$P.value, 
#          method = 'bonferroni')
# 
# 
# GO_cell_func = enriched_chrxxi_inversion$GO_Cellular_Component_2025%>% 
#   as_tibble() 
# p.adjust(GO_cell_func$P.value, 
#          method = 'bonferroni')
# 
# 
# GO_bio_func = enriched_chrxxi_inversion$GO_Biological_Process_2025%>% 
#   as_tibble() 
# p.adjust(GO_bio_func$P.value, 
#          method = 'bonferroni')
# 
# KEGG = enriched_chrxxi_inversion$KEGG_2021_Human%>% 
#   as_tibble() 
# p.adjust(KEGG$P.value, 
#          method = 'bonferroni')

## ******
PhenGen = enriched_chrxxi_inversion$PhenGenI_Association_2021%>% 
  as_tibble() 
PhenGen_pdj = p.adjust(PhenGen$P.value, 
         method = 'bonferroni')

PhenGen = bind_cols(PhenGen, 
                    PhenGen_pdj) %>% 
  rename(Bonferroni_adj_Pval = ...10)

PhenGen %>% 
  filter(Bonferroni_adj_Pval <= 0.05) %>% 
  View()

# 
# reactome = enriched_chrxxi_inversion$Reactome_Pathways_2024%>% 
#   as_tibble() 
# p.adjust(reactome$P.value, 
#          method = 'bonferroni')
# 
# 
# transfac = enriched_chrxxi_inversion$TRANSFAC_and_JASPAR_PWMs%>% 
#   as_tibble() 
# p.adjust(transfac$P.value, 
#          method = 'bonferroni')
# 
# 
# wiki_mouse = enriched_chrxxi_inversion$WikiPathways_2024_Mouse%>% 
#   as_tibble() 
# p.adjust(wiki_mouse$P.value, 
#          method = 'bonferroni')
# 
# wiki_human = enriched_chrxxi_inversion$WikiPathways_2024_Human%>% 
#   as_tibble() 
# p.adjust(wiki_human$P.value, 
#          method = 'bonferroni')



# GO biocarta pathways ----------------------------------------------------
#### *****************
dbs <- c("BioCarta_2013", 
         "BioCarta_2015",
         "BioCarta_2016",
         'BioPlanet_2019', 
         'BioPlex_2017')


biocarta_enriched <- enrichr(chrxxi_inversion_gene_names, 
                                     dbs, 
                                     background = stickle_background, 
                                     include_overlap = T)

biocarta_2013 = biocarta_enriched$BioCarta_2013 %>% 
  as_tibble() %>% 
  filter(Adjusted.P.value <= 0.05)


# GO diabetes perturb -----------------------------------------------------
#### *****************
dbs = c('Diabetes_Perturbations_GEO_2022')

diabetes_enriched <- enrichr(chrxxi_inversion_gene_names, 
                             dbs, 
                             background = stickle_background, 
                             include_overlap = T)

go_diabetes = diabetes_enriched$Diabetes_Perturbations_GEO_2022 %>% 
  as_tibble() %>% 
  filter(Adjusted.P.value <= 0.05)


# enrichr pathways --------------------------------------------------------
## ************88
dbs = c('Enrichr_Libraries_Most_Popular_Genes',
        'Enrichr_Submissions_TF-Gene_Coocurrence',
        'Enrichr_Users_Contributed_Lists_2020')

enrichr_enriched <- enrichr(chrxxi_inversion_gene_names, 
                             dbs, 
                             background = stickle_background, 
                             include_overlap = T)

enrichr_enriched$`Enrichr_Submissions_TF-Gene_Coocurrence` %>% 
  as_tibble() %>% 
  filter(Adjusted.P.value <= 0.05) 

# Go bio processes ----------------------------------------------------
## **************
dbs = c('GO_Biological_Process_2021',
        'GO_Biological_Process_2023',
        'GO_Biological_Process_2025')

GO_BIO_enriched <- enrichr(chrxxi_inversion_gene_names, 
                            dbs, 
                            background = stickle_background, 
                            include_overlap = T)

GO_BIO_enriched = GO_BIO_enriched$GO_Biological_Process_2021 %>% 
  as_tibble() %>% 
  filter(Adjusted.P.value <= 0.05) 


# GO cell processes -------------------------------------------------------
dbs = c('GO_Cellular_Component_2021',
        'GO_Cellular_Component_2023',
        'GO_Cellular_Component_2025')

GO_CELL_enriched <- enrichr(chrxxi_inversion_gene_names, 
                           dbs, 
                           background = stickle_background, 
                           include_overlap = T)

GO_CELL_enriched1 = GO_CELL_enriched$GO_Cellular_Component_2021 %>% 
  as_tibble() %>% 
  filter(Adjusted.P.value <= 0.05) 
GO_CELL_enriched2 = GO_CELL_enriched$GO_Cellular_Component_2023 %>% 
  as_tibble() %>% 
  filter(Adjusted.P.value <= 0.05) 
GO_CELL_enriched3 = GO_CELL_enriched$GO_Cellular_Component_2025 %>% 
  as_tibble() %>% 
  filter(Adjusted.P.value <= 0.05) 

# GO molecular function ---------------------------------------------------
dbs = c('GO_Molecular_Function_2021',
        'GO_Molecular_Function_2023',
        'GO_Molecular_Function_2025')

GO_MOL_enriched <- enrichr(chrxxi_inversion_gene_names, 
                            dbs, 
                            background = stickle_background, 
                            include_overlap = T)

GO_MOL_enriched1 = GO_MOL_enriched$GO_Molecular_Function_2021 %>% 
  as_tibble() %>% 
  filter(Adjusted.P.value <= 0.05) 
GO_MOL_enriched2 = GO_MOL_enriched$GO_Molecular_Function_2023 %>% 
  as_tibble() %>% 
  filter(Adjusted.P.value <= 0.05) 
GO_MOL_enriched3 = GO_MOL_enriched$GO_Molecular_Function_2025 %>% 
  as_tibble() %>% 
  filter(Adjusted.P.value <= 0.05) 


# GO GTeX database --------------------------------------------------------
dbs = c('GTEx_Aging_Signatures_2021',
        'GTEx_Tissues_V8_2023')

GO_GTEx_enriched <- enrichr(chrxxi_inversion_gene_names, 
                           dbs, 
                           background = stickle_background, 
                           include_overlap = T)

GO_GTEx_enriched$GTEx_Aging_Signatures_2021 %>% 
  as_tibble() %>% 
  filter(Adjusted.P.value <= 0.05) 
GO_GTEx_enriched$GTEx_Tissues_V8_2023 %>% 
  as_tibble() %>% 
  filter(Adjusted.P.value <= 0.05) 


# GO GWAS CATALOG ---------------------------------------------------------
dbs = c('GWAS_Catalog_2019',
        'GWAS_Catalog_2023', 
        'GWAS_Catalog_2025')

GO_GWAS_enriched <- enrichr(chrxxi_inversion_gene_names, 
                            dbs, 
                            background = stickle_background, 
                            include_overlap = T)

GO_GWAS1 = GO_GWAS_enriched$GWAS_Catalog_2019 %>% 
  as_tibble() %>% 
  filter(Adjusted.P.value <= 0.05) 
# GO_GWAS_enriched$GWAS_Catalog_2023 %>% 
#   as_tibble() %>% 
#   filter(Adjusted.P.value <= 0.05) 
GO_GWAS2 = GO_GWAS_enriched$GWAS_Catalog_2025 %>% 
  as_tibble() %>% 
  filter(Adjusted.P.value <= 0.05) 


# KEA database GO ---------------------------------------------------------


dbs = c('KEA_2013',
        'KEA_2015')

GO_KEA_enriched <- enrichr(chrxxi_inversion_gene_names, 
                            dbs, 
                            background = stickle_background, 
                            include_overlap = T)

GO_KEA_enriched$KEA_2013 %>% 
  as_tibble() %>% 
  filter(Adjusted.P.value <= 0.05) 
GO_KEA_enriched$KEA_2015 %>% 
  as_tibble() %>% 
  filter(Adjusted.P.value <= 0.05) 


# KEGG database GO --------------------------------------------------------
## ************
dbs = c('KEGG_2013',
        'KEGG_2015', 
        'KEGG_2016')

GO_KEGG_enriched <- enrichr(chrxxi_inversion_gene_names, 
                           dbs, 
                           background = stickle_background, 
                           include_overlap = T)

GO_KEGG1 = GO_KEGG_enriched$KEGG_2013 %>% 
  as_tibble() %>% 
  filter(Adjusted.P.value <= 0.05) 
# GO_KEGG_enriched$KEGG_2015 %>% 
#   as_tibble() %>% 
#   filter(Adjusted.P.value <= 0.05) 
# GO_KEGG_enriched$KEGG_2016 %>% 
#   as_tibble() %>% 
#   filter(Adjusted.P.value <= 0.05) 


# GO panther database -----------------------------------------------------
dbs = c('Panther_2015',
        'Panther_2016', 
        'PerturbAtlas')

GO_Panther_enriched <- enrichr(chrxxi_inversion_gene_names, 
                            dbs, 
                            background = stickle_background, 
                            include_overlap = T)

GO_Panther_enriched$Panther_2015 %>% 
  as_tibble() %>% 
  filter(Adjusted.P.value <= 0.05) 
GO_Panther_enriched$Panther_2016 %>% 
  as_tibble() %>% 
  filter(Adjusted.P.value <= 0.05) 

GO_Panther_enriched$PerturbAtlas %>% 
  as_tibble() %>% 
  filter(Adjusted.P.value <= 0.05) 


# Reactome GO -------------------------------------------------------------

dbs = c('Reactome_2022',
        'Reactome_Pathways_2024', 
        'The_Kinase_Library_2024', 
        'Transcription_Factor_PPIs')

GO_REACT_enriched <- enrichr(chrxxi_inversion_gene_names, 
                               dbs, 
                               background = stickle_background, 
                               include_overlap = T)

GO_REACT_enriched$Reactome_2022 %>% 
  as_tibble() %>% 
  filter(Adjusted.P.value <= 0.05) 
GO_REACT_enriched$Reactome_Pathways_2024 %>% 
  as_tibble() %>% 
  filter(Adjusted.P.value <= 0.05) 
GO_REACT_enriched$The_Kinase_Library_2024 %>% 
  as_tibble() %>% 
  filter(Adjusted.P.value <= 0.05) 
GO_REACT_enriched$Transcription_Factor_PPIs %>% 
  as_tibble() %>% 
  filter(Adjusted.P.value <= 0.05) 


# GO Wikipathways ---------------------------------------------------------
dbs = c('WikiPathways_2013',
        'WikiPathways_2015', 
        'WikiPathways_2016', 
        'WikiPathways_2019_Human', 
        'WikiPathways_2019_Mouse', 
        'WikiPathway_2021_Human', 
        'WikiPathway_2023_Human', 
        'WikiPathways_2024_Human', 
        'WikiPathways_2024_Mouse')

GO_WIKI_enriched <- enrichr(chrxxi_inversion_gene_names, 
                             dbs, 
                             background = stickle_background, 
                             include_overlap = T)

GO_WIKI_enriched$WikiPathways_2013 %>% 
  as_tibble() %>% 
  filter(Adjusted.P.value <= 0.05) 
GO_WIKI_enriched$WikiPathways_2015 %>% 
  as_tibble() %>% 
  filter(Adjusted.P.value <= 0.05) 
GO_WIKI_enriched$WikiPathways_2016 %>% 
  as_tibble() %>% 
  filter(Adjusted.P.value <= 0.05) 
GO_WIKI_enriched$WikiPathways_2019_Human %>% 
  as_tibble() %>% 
  filter(Adjusted.P.value <= 0.05) 
GO_WIKI_enriched$WikiPathways_2019_Mouse %>% 
  as_tibble() %>% 
  filter(Adjusted.P.value <= 0.05) 
GO_WIKI_enriched$WikiPathway_2021_Human %>% 
  as_tibble() %>% 
  filter(Adjusted.P.value <= 0.05) 
GO_WIKI_enriched$WikiPathway_2023_Human %>% 
  as_tibble() %>% 
  filter(Adjusted.P.value <= 0.05) 
GO_WIKI_enriched$WikiPathways_2024_Human %>% 
  as_tibble() %>% 
  filter(Adjusted.P.value <= 0.05) 
GO_WIKI_enriched$WikiPathways_2024_Mouse %>% 
  as_tibble() %>% 
  filter(Adjusted.P.value <= 0.05) 


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
          position) 
# %>% 
#   distinct(gene_name, 
#            .keep_all = T) %>% 
#   filter(!grepl('ENSG', 
#                 gene_name)) %>% 
#   filter(feature == 'gene')

ASHN_FST_out_genes %>% 
  write_csv('ASHN_FST_0.5%_outlier_genes_NEW.csv')

ASHN_FST_out_genes %>% 
  select(gene_name) %>% 
  write_tsv('ASHN_FST_0.5%_outlier_gene_names_only_NEW.tsv', 
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
          position) 
# %>% 
#   distinct(gene_name, 
#            .keep_all = T) %>% 
#   filter(!grepl('ENSG', 
#                 gene_name)) %>% 
#   filter(feature == 'gene')

MYV_FST_out_genes %>% 
  write_csv('MYV_FST_0.5%_outlier_genes_NEW.csv')

MYV_FST_out_genes %>% 
  select(gene_name) %>% 
  write_tsv('MYV_FST_0.5%_outlier_gene_names_only_NEW.tsv', 
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
          position) 
# %>% 
#   distinct(gene_name, 
#            .keep_all = T) %>% 
#   filter(!grepl('ENSG', 
#                 gene_name)) %>% 
  # filter(feature == 'gene')

SKR_FST_out_genes %>% 
  write_csv('SKR_FST_0.5%_outlier_genes_NEW.csv')

SKR_FST_out_genes %>% 
  select(gene_name) %>% 
  write_tsv('SKR_FST_0.5%_outlier_gene_names_only_NEW.tsv', 
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
          position) 

# %>% 
#   distinct(gene_name, 
#            .keep_all = T) %>% 
#   filter(!grepl('ENSG', 
#                 gene_name)) %>% 
#   filter(feature == 'gene')

GTS_CSWY_FST_out_genes %>% 
  write_csv('GTS_CSWY_FST_0.5%_outlier_genes_NEW.csv')


GTS_CSWY_FST_out_genes %>% 
  select(gene_name) %>% 
  write_tsv('GTS_CSWY_FST_0.5%_outlier_gene_names_only_NEW.tsv', 
            col_names = F)



# WC comparison  ----------------------------------------------------------

WC_div_snps = read_csv('~/Parsons_Postdoc/Stickleback_Genomic/vcf_filter/WC_FST_Outliers_WARM_COLD_COMBO.csv') %>% 
  # dplyr::select(-NMISS, 
  #        -FST) %>%
  rename(position = POS, 
         chromosome = CHR, 
         FST = FST) %>% 
  # filter(value == 'Outlier') %>% 
  arrange(chromosome, 
          position) %>% 
  group_by(chromosome)
## create a 100bp window around the snp identified
## settting the start and end range for the methylation data
## this number can be as big or small as you want depending
## on how liberal you want it

WC_div_window = WC_div_snps %>%
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

setDT(WC_div_window)
setDT(gene_annotation)

setkey(WC_div_window, 
       chromosome, 
       start, 
       end)

WC_gene_overlap = foverlaps(gene_annotation,
                                  WC_div_window,
                                  # by.x = start,
                                  # by.y = end,
                                  type="any")


WC_gene_overlap_tib = as_tibble(WC_gene_overlap) %>% 
  na.omit() %>% 
  filter(chromosome != 'chrUn') %>% 
  arrange(chromosome, 
          position)

gene_name_1 = WC_gene_overlap_tib %>% 
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

gene_name_2 = WC_gene_overlap_tib %>% 
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


WC_FST_out_genes = bind_rows(gene_name_1, 
                                   gene_name_2) %>% 
  arrange(chromosome, 
          position) 
# %>% 
  # distinct(gene_name, 
  #          .keep_all = T) %>% 
  # filter(!grepl('ENSG', 
  #               gene_name)) 
# %>% 
#   filter(feature == 'gene')

WC_FST_out_genes %>% 
  write_csv('WC_FST_0.5%_outlier_genes.csv')


WC_FST_out_genes %>% 
  select(gene_name) %>% 
  write_tsv('WC_FST_0.5%_outlier_gene_names_only.tsv', 
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
## create a 100bp window ar

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
# 
# as_tibble(ASHN_gene_overlap) %>% 
#   na.omit() %>% 
#   filter(chromosome != 'chrUn') %>% 
#   arrange(chromosome, 
#           position) %>% 
#   filter(feature == 'gene')

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
    filter(!grepl('ENSG',
                  gene_name)) %>% 
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

as_tibble(MYV_gene_overlap) %>% 
  na.omit() %>% 
  filter(chromosome != 'chrUn') %>% 
  arrange(chromosome, 
          position) %>% 
  filter(feature == 'gene')

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
          position) 
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




# All loci gene annotation ------------------------------------------------

ASHN_Fst_clean = read_csv('ASHN_TOP_DAWG_Fst_clean.csv')
MYV_Fst_clean = read_csv('MYV_TOP_DAWG_Fst_clean.csv') 
SKR_Fst_clean = read_csv('SKR_TOP_DAWG_Fst_clean.csv') 
GTS_CSWY_Fst_clean = read_csv('GTS_CSWY_TOP_DAWG_Fst_clean.csv') 



# ASHN all loci genes -----------------------------------------------------

ASHN_snps = read_csv('ASHN_TOP_DAWG_Fst_clean.csv') %>% 
  dplyr::select(-NMISS,
         -FST) %>%
  rename(position = POS, 
         chromosome = CHR, 
         FST = FST_zero) %>% 
  # filter(value == 'Outlier') %>% 
  arrange(chromosome, 
          position) %>% 
  group_by(chromosome)
## create a 100bp window around the snp identified
## settting the start and end range for the methylation data
## this number can be as big or small as you want depending
## on how liberal you want it

ASHN_window = ASHN_snps %>%
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

setDT(ASHN_window)
setDT(gene_annotation)

setkey(ASHN_window, 
       chromosome, 
       start, 
       end)

ASHN_gene_overlap = foverlaps(gene_annotation,
                            ASHN_window,
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


ASHN_FST_genes = bind_rows(gene_name_1, 
                             gene_name_2) %>% 
  arrange(chromosome, 
          position) 
# %>% 
# distinct(gene_name, 
#          .keep_all = T) %>% 
# filter(!grepl('ENSG', 
#               gene_name)) 
# %>% 
#   filter(feature == 'gene')


ASHN_FST_genes %>% 
  filter(gene_name == 'braf')

ASHN_FST_genes %>% 
  filter(gene_name == 'prkcg')

ASHN_FST_genes %>% 
  write_csv('ASHN_FST_genes.csv')


ASHN_FST_genes %>% 
  select(gene_name) %>% 
  write_tsv('ASHN_FST_gene_names_only.tsv', 
            col_names = F)


# MYV all loci  -----------------------------------------------------------
MYV_snps = read_csv('MYV_TOP_DAWG_Fst_clean.csv') %>% 
  dplyr::select(-NMISS,
                -FST) %>%
  rename(position = POS, 
         chromosome = CHR, 
         FST = FST_zero) %>% 
  # filter(value == 'Outlier') %>% 
  arrange(chromosome, 
          position) %>% 
  group_by(chromosome)
## create a 100bp window around the snp identified
## settting the start and end range for the methylation data
## this number can be as big or small as you want depending
## on how liberal you want it

MYV_window = MYV_snps %>%
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

setDT(MYV_window)
setDT(gene_annotation)

setkey(MYV_window, 
       chromosome, 
       start, 
       end)

MYV_gene_overlap = foverlaps(gene_annotation,
                              MYV_window,
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


MYV_FST_genes = bind_rows(gene_name_1, 
                           gene_name_2) %>% 
  arrange(chromosome, 
          position) 
# %>% 
# distinct(gene_name, 
#          .keep_all = T) %>% 
# filter(!grepl('ENSG', 
#               gene_name)) 
# %>% 
#   filter(feature == 'gene')


MYV_FST_genes %>% 
  filter(gene_name == 'braf')

MYV_FST_genes %>% 
  filter(gene_name == 'prkcg')

MYV_FST_genes %>% 
  write_csv('MYV_FST_genes.csv')


MYV_FST_genes %>% 
  select(gene_name) %>% 
  write_tsv('MYV_FST_gene_names_only.tsv', 
            col_names = F)


# SKR ---------------------------------------------------------------------

SKR_snps = read_csv('SKR_TOP_DAWG_Fst_clean.csv') %>% 
  dplyr::select(-NMISS,
                -FST) %>%
  rename(position = POS, 
         chromosome = CHR, 
         FST = FST_zero) %>% 
  # filter(value == 'Outlier') %>% 
  arrange(chromosome, 
          position) %>% 
  group_by(chromosome)
## create a 100bp window around the snp identified
## settting the start and end range for the methylation data
## this number can be as big or small as you want depending
## on how liberal you want it

SKR_window = SKR_snps %>%
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

setDT(SKR_window)
setDT(gene_annotation)

setkey(SKR_window, 
       chromosome, 
       start, 
       end)

SKR_gene_overlap = foverlaps(gene_annotation,
                              SKR_window,
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


SKR_FST_genes = bind_rows(gene_name_1, 
                           gene_name_2) %>% 
  arrange(chromosome, 
          position) 
# %>% 
# distinct(gene_name, 
#          .keep_all = T) %>% 
# filter(!grepl('ENSG', 
#               gene_name)) 
# %>% 
#   filter(feature == 'gene')


SKR_FST_genes %>% 
  filter(gene_name == 'braf')

SKR_FST_genes %>% 
  filter(gene_name == 'prkcg')

SKR_FST_genes %>% 
  write_csv('SKR_FST_genes.csv')


SKR_FST_genes %>% 
  select(gene_name) %>% 
  write_tsv('SKR_FST_gene_names_only.tsv', 
            col_names = F)


# gts_cswy alll loci and genes --------------------------------------------

GTS_CSWY_snps = read_csv('GTS_CSWY_TOP_DAWG_Fst_clean.csv') %>% 
  dplyr::select(-NMISS,
                -FST) %>%
  rename(position = POS, 
         chromosome = CHR, 
         FST = FST_zero) %>% 
  # filter(value == 'Outlier') %>% 
  arrange(chromosome, 
          position) %>% 
  group_by(chromosome)
## create a 100bp window around the snp identified
## settting the start and end range for the methylation data
## this number can be as big or small as you want depending
## on how liberal you want it

GTS_CSWY_window = GTS_CSWY_snps %>%
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

setDT(GTS_CSWY_window)
setDT(gene_annotation)

setkey(GTS_CSWY_window, 
       chromosome, 
       start, 
       end)

GTS_CSWY_gene_overlap = foverlaps(gene_annotation,
                              GTS_CSWY_window,
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


GTS_CSWY_FST_genes = bind_rows(gene_name_1, 
                           gene_name_2) %>% 
  arrange(chromosome, 
          position) 
# %>% 
# distinct(gene_name, 
#          .keep_all = T) %>% 
# filter(!grepl('ENSG', 
#               gene_name)) 
# %>% 
#   filter(feature == 'gene')


GTS_CSWY_FST_genes %>% 
  filter(gene_name == 'braf')

GTS_CSWY_FST_genes %>% 
  filter(gene_name == 'prkcg')

GTS_CSWY_FST_genes %>% 
  write_csv('GTS_CSWY_FST_genes.csv')


GTS_CSWY_FST_genes %>% 
  select(gene_name) %>% 
  write_tsv('GTS_CSWY_FST_gene_names_only.tsv', 
            col_names = F)






