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

map_data = read_tsv('stickleback_maf0.05_ldpruned_filtered.map', 
                    col_names = c('CHR', 
                                  'SNP', 
                                  'GPOS', 
                                  'POS'))

ped_data = read_table('stickleback_maf0.05_ldprunded_filtered.raw',
                      col_names = c('PopulationID', 
                                    'IndividualID', 
                                    'MaternalID', 
                                    'PaternalID', 
                                    'Sex', 
                                    'Phenotype', 
                                    map_data$SNP))

tped = Create_tped(ped = ped_data, 
            map = map_data)

# tped %>%
#   write_tsv('stickleback_maf0.05_ldpruned_filtered.tped', 
#             col_names = F)
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



# Lostruct on full chr XXI -------------------------------------------------


## Need to plot the MDS structure for the regions 
## flanking the potential inversion

theme_set(theme_bw())

Normal_data = combo_data %>% 
  filter(!window %in% c('83',
                        '84', 
                        '85', 
                        '86', 
                        '87', 
                        '88', 
                        '89', 
                        '90', 
                        '91', 
                        '92',
                        '93')) 
label = rep('Normal_region', 
     nrow(Normal_data))

Normal_data = cbind(label, 
                    Normal_data) %>% 
  as_tibble() %>% 
  rename(group = 1)%>% 
  dplyr::select(group, 
                X4, 
                window, 
                MDS_Points1, 
                MDS_Points2)
# afvaper_region = combo_data %>%
#   filter(window >= '81',
#          window <= '95')
# MDS_points_windows = afvaper_region %>% 
#   dplyr::select(X4, 
#                 window, 
#                 MDS_Points1, 
#                 MDS_Points2)
# 
# MDS_points_windows %>% 
#   write_csv('AFvaper_region_lostruct_survey.csv')

afvaper_region = read_csv('AFvaper_region_lostruct_survey.csv') %>% 
  filter(group == 'Afvaper_region')

labeled_data = bind_rows(Normal_data, 
                         afvaper_region)

# window_distance = afvaper_region %>%
#   dplyr::select(contains('V'))

labeled_data %>% 
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
MDS_outliers = Outlier_hunter(data = labeled_data, 
                              sd_percentile = 3) 

Outlier_plots(outlier_data = MDS_outliers, 
              normal_data = Normal_data)



# Calculate window size - 50 snp windows ----------------------------------

Normal_data_win_size = combo_data %>% 
  dplyr::select(X1, 
                window, 
                mean_window) %>% 
  group_by(window) %>% 
  mutate(mean_window_Mb = mean_window/1000000)

Normal_data_win_size$mean_window_Mb

outlier_data_win_size = combo_data %>% 
  filter(window %in% c('83',
                        '84', 
                        '85', 
                        '86', 
                        '87', 
                        '88', 
                        '89', 
                        '90', 
                        '91', 
                        '92')) %>% 
  dplyr::select(X1, 
                window, 
                mean_window) %>% 
  group_by(window) %>% 
  mutate(mean_window_Mb = mean_window/1000000)


##s
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

stickle_plot = bind_cols(identifiers, 
                         stickle_plot)

stickle_plot = mutate(.data = stickle_plot,
                      Location = as.factor(case_when(
                        population == 'ASHNC' ~ 'Áshildarholtsvatn',
                        population == 'ASHNW' ~ 'Áshildarholtsvatn',
                        population == 'CSWY' ~ 'Garðsvatn',
                        population == 'GTS' ~ 'Grettislaug',
                        population == 'MYVC' ~ 'Mývatn',
                        population == 'MYVW' ~ 'Mývatn',
                        population == 'SKRC' ~ 'Sauðárkrókur',
                        population == 'SKRW' ~ 'Sauðárkrókur')))

stickle_plot = mutate(.data = stickle_plot, 
                      Type = as.factor(case_when(
                        population == 'ASHNC' ~ 'Ambient',
                        population == 'ASHNW' ~ 'Geothermal',
                        population == 'CSWY' ~ 'Ambient',
                        population == 'GTS' ~ 'Geothermal',
                        population == 'MYVC' ~ 'Ambient',
                        population == 'MYVW' ~ 'Geothermal',
                        population == 'SKRC' ~ 'Ambient',
                        population == 'SKRW' ~ 'Geothermal'
                      )))


theme_set(theme_bw())

location_cols = c('#06d6a0',
                  '#264653',
                  '#219ebc',
                  '#d62828',
                  '#5f0f40')


stickleback_pca = stickle_plot %>%
  # arrange(population) %>% 
  # group_by(population) %>% 
  ggplot(aes(x = PC1, 
             y = PC2))+
  geom_point(aes(col = Location, 
                 shape = Type),
             size = 3)+
  # geom_point(aes(col = population),
  #            size = 2)+
  # scale_color_manual(values = cold_warm_cols)+
  scale_color_manual(values = location_cols)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.title = element_text(size = 15), 
        axis.text = element_text(size = 15), 
        axis.ticks = element_line(size = 1), 
        plot.title = element_text(size = 15, 
                                  hjust = 0), 
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14)) +
  labs(x = 'Principal component 1 (28.96%)',
       y = 'Principal component 2 (10.65%)', 
       col = 'Populations')

stickleback_pca


ggsave(file = 'PCA_Afvaper_Inversion_region_lostruct_outlier.tiff', 
       path = '~/Parsons_Postdoc/Stickleback_Genomic/Figures/', 
       plot = stickleback_pca, 
       dpi = 'retina', 
       units = 'cm', 
       width = 20.0, 
       height = 13)






# DAPCA chr21 inversion region --------------------------------------------

library(adegenet)

map = read_tsv('chr21_inversion_region.map', 
               col_names = F) %>% 
  rename(CHR = X1, 
         SNP = X2, 
         GPOS = X3, 
         POS = X4)

chr21_inversion_data = read_tsv('chr21_inversion_region.ped', 
                    col_names = c('#FamilyID', 
                                  'IndividualID', 
                                  'PaternalID', 
                                  'MaternalID', 
                                  'Sex', 
                                  'Phenotype', 
                                  map$SNP)) %>% 
  dplyr::select(`#FamilyID`, 
                IndividualID, 
                contains('chr_'))


## Need to change to a data frame from a tibble
## otherwise we'll get a weird dimension names error
chr21_inversion_data = as.data.frame(chr21_inversion_data)

genind = df2genind(chr21_inversion_data[,3:length(chr21_inversion_data)], 
                       ploidy = 2, 
                       ind.names = chr21_inversion_data[,2], 
                       sep = '\t', 
                       pop = chr21_inversion_data[,1])
set.seed(666)
chr21_inversion_clust = find.clusters(genind, 
                                      max.n.clust = 10)

## Figured it out....i'm a coding GOD
dapc_table = table(pop(genind), 
                   chr21_inversion_clust$grp)

## write the table to see if the strong clustering 
## relates chr
# as.data.frame(dapc_table) %>% 
#   write_tsv('DAPC_output_three_clusters.txt')

table.value(table(pop(genind), 
                  chr21_inversion_clust$grp), 
            col.lab=paste("inf", 1:4),
            row.lab=paste("ori", 1:37))

chr21_inversion_clust$grp %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  as_tibble() %>% 
  rename(IndividualID = 1, 
         dapc_group_assignment = 2) %>% 
  write_csv('DAPCA_analysis_chr21_inversion_region.csv')


## be careful and pay attention
## the number of PC axes you use with this part of the analysis
## will largely change the dapc output
dapc_chr21_region = dapc(genind, 
                chr21_inversion_clust$grp)

## contribution of individuals to the 2 LDs

as.data.frame(dapc_chr21_region$ind.coord) %>%
  as_tibble() %>%
  write_csv('Dapc_individual_LD_coordinates.csv')

## contribution of the groups to the 3 LDs

as.data.frame(dapc_chr21_region$grp.coord) %>%
  as_tibble() %>%
  write_csv('Dapc_group_LD_coordinates.csv')

## contribution of the individual loci to the 3 LDs

as.data.frame(dapc_chr21_region$var.contr) %>%
  as_tibble() %>%
  write_csv('Dapc_SNP_LD_coordinates.csv')

## eigenvalues
dapc_chr21_region$eig

## individuals are dots
## groups are inertia ellipses
scatter(dapc_chr21_region)


pal = c('#2E4159',
        '#4968A6',
        # '#F29E38',
        '#F23E2E')

scatter(dapc_chr21_region, 
        scree.da=FALSE, 
        bg="white", 
        pch=20, 
        cell=0, 
        cstar=0, 
        col=pal,
        solid=.4,
        cex=3,
        clab=0, 
        leg=FALSE)
# txt.leg=paste("Cluster",1:3))

scatter(dapc_chr21_region,
        1,
        1, 
        col=pal, 
        bg="white",
        scree.da=FALSE, 
        legend=F, 
        solid=.4)

# Dapc ggplot mashup ------------------------------------------------------

# individual_coords = read_csv('~/Charr_Adaptive_Introgression/Charr_Project_1/GeneticData/Lab_dapc_individual_LD_coordinates.csv')
group_coords = read_csv('~/Charr_Adaptive_Introgression/Charr_Project_1/GeneticData/Lab_dapc_group_LD_coordinates.csv')
snp_coords = read_csv('~/Charr_Adaptive_Introgression/Charr_Project_1/GeneticData/Lab_dapc_SNP_LD_coordinates.csv')
individual_coords = read_csv('~/Charr_Adaptive_Introgression/Charr_Project_1/GeneticData/Lab_dapc_individual_coords.csv')
individual_groups = read_csv('~/Charr_Adaptive_Introgression/Charr_Project_1/GeneticData/dapc_Individual_group_assignment.csv')
# indiv = read_csv('~/Charr_Adaptive_Introgression/Charr_Project_1/GeneticData/Lab_dapc_individual_LD_coordinates.csv')
# lab_metadata = lab_data %>%
#   as_tibble() %>%
#   dplyr::select(`#FamilyID`,
#                 IndividualID) %>%
#   rename(Population = `#FamilyID`)
# coords = read_csv('~/Charr_Adaptive_Introgression/Charr_Project_1/SampleSiteData/LabSampleSites_Coord_1June2020.csv')
# 
# lab_metadata = full_join(lab_metadata,
#           coords,
#           by = 'Population')
# 
# bind_cols(lab_metadata,
#                      individual_coords) %>%
#   write_csv('~/Charr_Adaptive_Introgression/Charr_Project_1/GeneticData/Lab_dapc_individual_coords.csv')

individual_data = bind_cols(individual_coords, 
                            individual_groups)

individual_data$dapc_group_assignment = as.character(individual_data$dapc_group_assignment)

library(ggforce)
library(patchwork)

pal = c('#2E4159',
        '#4968A6',
        '#F23E2E')

theme_set(theme_bw())

dapc_plot1 = ggplot(data = individual_data,
                    aes(x = LD1, 
                        y = LD2))+
  # geom_point(aes(col = dapc_group_assignment), 
  #            size = 1.5)+
  geom_point(aes(col = Lat...4),
             size = 1.5)+
  # geom_mark_ellipse(expand = 0, 
  #                   size = 2,
  #                  aes(col = dapc_group_assignment))+
  geom_mark_ellipse(expand = 0, 
                    size = 1, 
                    aes(group = dapc_group_assignment))+
  # scale_color_manual(values = pal)+
  labs(x = 'Linear discriminant axis 1', 
       y = 'Linear discriminant axis 2', 
       color = 'Latitude')+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.title.y = element_blank(),
        axis.text =  element_text(size = 12),
        axis.ticks = element_line(size = 1),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))
dapc_plot1

dapc_plot2 = ggplot(data = individual_data,
                    aes(x = LD1, 
                        y = LD2))+
  geom_point(aes(col = dapc_group_assignment),
             size = 1.5)+
  geom_mark_ellipse(expand = 0,
                    size = 1,
                    aes(group = dapc_group_assignment))+
  scale_color_manual(values = pal)+
  labs(x = 'Linear discriminant axis 1', 
       y = 'Linear discriminant axis 2', 
       color = 'Latitude')+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text =  element_text(size = 12),
        axis.ticks = element_line(size = 1),
        legend.position = 'none')

dapc_plot3 = ggplot(data = individual_data, 
                    aes(x = LD1))+
  geom_histogram(col = 'black',
                 aes(fill = dapc_group_assignment))+
  scale_fill_manual(values = pal)+
  labs(x = 'Linear discriminant axis 1', 
       y = 'Count')+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text =  element_text(size = 12),
        axis.ticks = element_line(size = 1),
        legend.position = 'none')

dapc_plot4 = ggplot(data = individual_data, 
                    aes(x = LD2))+
  geom_histogram(col = 'black',
                 aes(fill = dapc_group_assignment))+
  scale_fill_manual(values = pal)+
  labs(x = 'Linear discriminant axis 2', 
       y = 'Count', 
       fill = 'Cluster')+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.title.y = element_blank(),
        axis.text =  element_text(size = 12),
        axis.ticks = element_line(size = 1), 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 12))

combo_plot = (dapc_plot2 + dapc_plot1)/(dapc_plot3 + dapc_plot4)

