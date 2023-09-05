##############################
## lostruct stickleback chr21
##
## Matt Brachmann (PhDMattyB)
##
## 30.08.2023
##
##############################

setwd('~/Parsons_Postdoc/Stickleback_Genomic/vcf_filter/')

library(tidyverse)
library(data.table)

ped_data = read_table('stickleback_maf0.05_ldprunded_filtered.raw',
                      col_names = c('PopulationID', 
                                    'IndividualID', 
                                    'MaternalID', 
                                    'PaternalID', 
                                    'Sex', 
                                    'Phenotype', 
                                    map_data$SNP))

map_data = read_tsv('stickleback_maf0.05_ldpruned_filtered.map', 
                    col_names = c('CHR', 
                                  'SNP', 
                                  'GPOS', 
                                  'POS'))

tped = Create_tped(ped = ped_data, 
            map = map_data)

tped %>%
  write_tsv('stickleback_maf0.05_ldpruned_filtered.tped', 
            col_names = F)
# 
snps = read_table('stickleback_maf0.05_ldpruned_filtered.tped', 
                  col_names = F)

df = snps %>% 
  dplyr::select(1, 
                5:length(snps)) %>% 
  filter(X1 == 'chr_XXI') %>% 
  dplyr::select(-X1) %>% 
  dplyr::select(-X5)

eigen = eigen_windows(df, 
                      win = 50, 
                      k = 5)
windist = pc_dist(eigen, 
                  npc = 5) %>% 
  as_tibble()

window_data = snps %>% 
  select(1:4) %>% 
  filter(X1 == 'chr_XXI') %>% 
  mutate(window = ceiling(row_number()/50)) %>% 
  group_by(window) %>% 
  mutate(mean_window = mean(X4)) %>% 
  distinct(mean_window, 
           .keep_all = T) %>% 
  filter(window %in% 1:nrow(windist))

combo_data = bind_cols(window_data, 
                       windist)


MDS_data = cmdscale(combo_data[7:length(combo_data)], 
                    eig = TRUE, 
                    k = 5)

MDS_points = MDS_data$points %>% 
  as_tibble() %>% 
  rename(MDS_Points1 = V1, 
         MDS_Points2 = V2, 
         MDS_points3 = V3, 
         MDS_points4 = V4, 
         MDS_points5 = V5)

combo_data = bind_cols(combo_data, 
                       MDS_points) %>% 
  group_by(window)

View(combo_data)
# plot the MDS structure  -------------------------------------------------


## Need to plot the MDS structure for the regions 
## flanking the potential inversion

## windows 83 and 93 are the start and end positions of the afvaper region
theme_set(theme_bw())

afvaper_region = combo_data %>%
  filter(window >= '81',
         window <= '95')
# MDS_points_windows = afvaper_region %>% 
#   dplyr::select(X4, 
#                 window, 
#                 MDS_Points1, 
#                 MDS_Points2)
# 
# MDS_points_windows %>% 
#   write_csv('AFvaper_region_lostruct_survey.csv')

MDS_points_windows = read_csv('AFvaper_region_lostruct_survey.csv')

window_distance = afvaper_region %>%
  dplyr::select(contains('V'))

MDS_points_windows %>% 
  ggplot(aes(x = MDS_Points1, 
             y = MDS_Points2,
             col = group))+
  # ggplot(aes(x = MDS_Points1, 
  #            y = MDS_Points2,
  #            col = rainbow(nrow(window_distance))))+
  geom_point(size = 3) +
  labs(x = 'MDS coordinate 1', 
       y = 'MDS coordinate 2')
# +
#   theme(legend.position = 'none')



## HOLY SHIT!! window 83-92 are MDS outliers!!!
MDS_outliers = Outlier_hunter(data = MDS_points_windows, 
               sd_percentile = 3) 

Normal_data = MDS_points_windows %>% 
  filter(window %in% c('81', 
                       '82', 
                       '93', 
                       '94', 
                       '95'))

Outlier_plots(outlier_data = MDS_outliers, 
              normal_data = Normal_data)



# PCA of the CHR21 region -------------------------------------------------

## Holy fuck just filter in plink with the --from --to flags to filter
## by the specific snp range
## the .raw file from RecodeA data can't be loaded by plink. 
ped_data = read_table('stickleback_maf0.05_ldpruned_filtered.ped',
                      col_names = c('PopulationID',
                                    'IndividualID',
                                    'MaternalID',
                                    'PaternalID',
                                    'Sex',
                                    'Phenotype',
                                    map_data$SNP))

chr21_ped_data = ped_data %>%
  # slice(-1) %>%
  dplyr::select(PopulationID,
                IndividualID,
                MaternalID,
                PaternalID,
                Sex,
                Phenotype,
                'chr_XXI_9963830':'chr_XXI_11370710')

chr21_map_data = read_tsv('stickleback_maf0.05_ldpruned_filtered.map',
                    col_names = c('CHR',
                                  'SNP',
                                  'GPOS',
                                  'POS')) %>%
  filter(CHR == 'chr_XXI') %>%
  filter(POS >= 9963830,
         POS <= 11370710)

chr21_ped_data %>%
  write_tsv('chr21_inversion_region.ped',
            col_names = F)
chr21_map_data %>%
  write_tsv('chr21_inversion_region.map',
            col_names = F)



# PCA on the putative inversion -------------------------------------------

library(pcadapt)


chr21_inversion = read.pcadapt('chr21_inversion_region.bed', 
                               type = 'bed')

chr21_pca = pcadapt::pcadapt(chr21_inversion, 
                               K = 10, 
                               method = 'mahalanobis', 
                               min.maf = 0.01)

plot(chr21_pca, 
     option = 'screeplot', 
     K = 10)


## percent variation explained by the three significant axes
chr21_pca$singular.values
sum(chr21_pca$singular.values)
(0.6293236/2.172823)*100
(0.2314772/2.172823)*100
## singular values should have already been square root transformed

chr21_pca_scores = as_tibble(chr21_pca$scores) %>%
  rename(PC1 = 1,
         PC2 = 2) %>%
  dplyr::select(PC1,
                PC2) %>%
  write_csv('chr21_PCA_scores.csv')

stickle_plot = read_csv('chr21_PCA_scores.csv')

identifiers = read_csv('stickleback_identifiers.csv')

