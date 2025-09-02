##############################
## stickleback genomics explore
##
## Matt Brachmann (PhDMattyB)
##
## 2022-11-23
##
##############################

library(pcadapt)
library(viridis)
# library(LEA)
library(reshape2)
library(devtools)
# install_github("jdstorey/qvalue")
library(qvalue)
library(tidyverse)
library(umap)
library(vcfR)
library(Rcpp)
# install_github('tavareshugo/windowscanr')
library(windowscanr)

# 
theme_set(theme_bw())

# setwd('C:/Stickleback_Genomic/vcf_filter/')
setwd('~/Parsons_Postdoc/Stickleback_Genomic/vcf_filter/')

##



# Tagging Exp 1 -----------------------------------------------------------
setwd('~/Parsons_Postdoc/Experiment1/')

read_csv('Experiment1_Elastomer_Tagging.csv') %>% 
  group_by(Population,
           tagid) %>% 
  distinct(.keep_all = T)


# Cross_Numbers -----------------------------------------------------------

setwd('~/Parsons_Postdoc/Experiment1')

cross_num = read_csv('Fish_Cross_Numbers_30.11.22.csv')

cross_num %>% 
  group_by(`P (m/f)`) %>% 
  filter(`P (m/f)` %in% c('acac', 
                          'awaw', 
                          'mcmc', 
                          'mwmw')) %>% 
  summarize(totalw = sum(totalw), 
            totalc = sum(totalc))

cross_num %>% 
  # filter(`P (m/f)`  == 'mcmc') %>% 
  filter(`P (m/f)` == 'awaw') %>% 
  dplyr::select(1:2, 
                totalw, 
                totalc) %>% 
  group_by(Cross) %>% 
  summarise(totalw = sum(totalw), 
            totalc = sum(totalc)) %>% 
  View()



# admixture inversion region ----------------------------------------------

inversion_map = read_tsv('chr21_inversion_region.map', 
                         col_names = c('Chromosome', 
                                       'SNP',
                                       'GPos', 
                                       'PPos')) %>% 
  Chr_convert() %>% 
  select(-Chromosome) %>% 
  select(chr_num, 
         SNP, 
         GPos, 
         PPos) %>% 
  write_tsv('chr21_inversion_map_fixed.map', 
            col_names = F)

inversion_ped = read_table('chr21_inversion_region.ped', 
                           col_names = F)



# admixture K4 chr21 inversion regions ---------------------------------------
## the cv error for K4 is the best. 
## check what k5 looks like as well. 

K4_Qval = read_table('chr21_inversion_region.4.Q', 
                     col_names = c('Q1', 
                                   'Q2', 
                                   'Q3', 
                                   'Q4'))

identifiers = read_csv('stickleback_identifiers.csv')
identifiers = identifiers %>% 
  rowid_to_column() %>% 
  rename(order = rowid)

K4_data = bind_cols(identifiers, 
                    K4_Qval)%>% 
  write_csv('Admixture_inversion_region_K4.csv')


K4_data = read_csv('Admixture_inversion_region_K4.csv')

K4_melted_data = melt(K4_data, 
                      id.vars = c('Plot_order',
                        'order',
                                  'population', 
                                  'individual_id')) %>% 
  as_tibble() 

K4_cols = c('#00798c',
            '#003d5b',
            '#edae49',
            '#d1495b')

k4_inversion_plot = ggplot(data = K4_melted_data, 
                           aes(x = reorder(Plot_order, 
                                           individual_id),
                               y = value, 
                               fill = variable), 
                           col = 'black')+
  geom_bar(stat = "identity", 
           width = 1, 
           col = 'black')+
  scale_fill_manual(values = K4_cols)+
  # scale_fill_manual(values = magma(n = 4))+
  labs(x = 'Individuals', 
       y = 'Ancestry proportion')+
  theme(axis.text.y = element_text(color = 'black', 
                                   size = 12),
        axis.title.y = element_text(color = 'black', 
                                    size = 14),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        ## can add xaxis labels if needed
        # axis.text.x = element_text(angle = 90,
        #                            hjust = 1,
        #                            vjust = -0.09,
        #                            size = 8,
        #                            color = 'black'),
        legend.position = 'none')+
  scale_x_discrete(guide = guide_axis(n.dodge = 5))+
  scale_y_continuous(expand = c(0,0))

k4_inversion_plot

ggsave('Chr21_inversion_admixture_K4.tiff',
       plot = k4_inversion_plot, 
       dpi = 'retina', 
       units = 'cm', 
       width = 35, 
       height = 10)


# admixture k3 inversion region -------------------------------------------
K3_Qval = read_table('chr21_inversion_region.3.Q', 
                     col_names = c('Q1', 
                                   'Q2', 
                                   'Q3'))

identifiers = read_csv('stickleback_identifiers.csv')
identifiers = identifiers %>% 
  rowid_to_column() %>% 
  rename(order = rowid)

K3_data = bind_cols(identifiers, 
                    K3_Qval)%>% 
  write_csv('Admixture_inversion_region_K3.csv')


K3_data = read_csv('Admixture_inversion_region_K3.csv')

K3_melted_data = melt(K3_data, 
                      id.vars = c('Plot_order',
                                  'order',
                                  'population', 
                                  'individual_id')) %>% 
  as_tibble() 

location_cols = c('#edae49',
                  '#d1495b',
                  '#00798c')
                  

K3_inversion_plot = ggplot(data = K3_melted_data, 
                           aes(x = reorder(Plot_order, 
                                           individual_id),
                               y = value, 
                               fill = variable), 
                           col = 'black')+
  geom_bar(stat = "identity", 
           width = 1, 
           col = 'black')+
  scale_fill_manual(values = location_cols)+
  # scale_fill_manual(values = magma(n = 4))+
  labs(x = 'Individuals', 
       y = 'Ancestry proportion')+
  theme(axis.text.y = element_text(color = 'black', 
                                   size = 12),
        axis.title.y = element_text(color = 'black', 
                                    size = 14),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        ## can add xaxis labels if needed
        # axis.text.x = element_text(angle = 90,
        #                            hjust = 1,
        #                            vjust = -0.09,
        #                            size = 8,
        #                            color = 'black'),
        legend.position = 'none')+
  scale_x_discrete(guide = guide_axis(n.dodge = 5))+
  scale_y_continuous(expand = c(0,0))

K3_inversion_plot

ggsave('admixture_k3_GTS_MYV_other.tiff',
       plot = admixture_k3_plot, 
       dpi = 'retina', 
       units = 'cm', 
       width = 35, 
       height = 10)



# admixture results -------------------------------------------------------
setwd('Parsons_Postdoc/Stickleback_Genomic/vcf_filter/')

## K3 CV error is 0.47767
## K4 CV error is 0.47944
## Pretty close. K3 makes sense based off of the PCA

K3_Qval = read_table('stickleback_maf0.05_ldpruned_filter_chr_fix.3.Q', 
                   col_names = c('Q1', 
                                 'Q2', 
                                 'Q3'))
# K3_Pval = read_table('stickleback_maf0.05_ldpruned_filter_chr_fix.3.P',
#                    col_names = c('P1',
#                                  'P2',
#                                  'P3'))

identifiers = read_csv('stickleback_identifiers.csv')
identifiers = identifiers %>% 
  rowid_to_column() %>% 
  rename(order = rowid)

K3_data = bind_cols(identifiers, 
                    K3_Qval)%>% 
  write_csv('Admixture_K3_Data.csv')

K3_data = read_csv('Admixture_K3_Data.csv')

K3_data %>% 
  group_by(population) %>% 
  summarize(n = n())



K3_melted_data = melt(K3_data, 
                    id.vars = c('Plot_order', 
                                'order',
                                'population', 
                                'individual_id')) %>% 
  as_tibble() 

K3_cols = c('#edae49',
            '#d1495b',
            '#00798c')

admixture_k3_plot = ggplot(data = K3_melted_data, 
       aes(x = reorder(Plot_order, 
                       individual_id),
           y = value, 
           fill = variable), 
       col = 'black')+
  geom_bar(stat = "identity", 
           width = 1, 
           col = 'black')+
  scale_fill_manual(values = K3_cols)+
  # scale_fill_manual(values = magma(n = 4))+
  labs(x = 'Individuals', 
       y = 'Ancestry proportion')+
  theme(axis.text.y = element_text(color = 'black', 
                                   size = 12),
        axis.title.y = element_text(color = 'black', 
                                    size = 14),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        ## can add xaxis labels if needed
        # axis.text.x = element_text(angle = 90,
        #                            hjust = 1,
        #                            vjust = -0.09,
        #                            size = 8,
        #                            color = 'black'),
        legend.position = 'none')+
  scale_x_discrete(guide = guide_axis(n.dodge = 5))+
  scale_y_continuous(expand = c(0,0))

admixture_k3_plot

ggsave('admixture_k3_GTS_MYV_other.tiff',
       plot = admixture_k3_plot, 
       dpi = 'retina', 
       units = 'cm', 
       width = 35, 
       height = 10)



# K4_cols = c('#edae49',
#             '#d1495b',
#             '#00798c')
# 
# K4_Qval = read_table('stickleback_maf0.05_ldpruned_filter_chr_fix.4.Q', 
#                      col_names = c('Q1', 
#                                    'Q2', 
#                                    'Q3', 
#                                    'Q4'))
# # K4_Pval = read_tsv('stickleback_maf0.05_ldpruned_filter_chr_fix.3.P', 
# #                    col_names = c('P1', 
# #                                  'P2', 
# #                                  'P3'))
# 
# identifiers = read_csv('stickleback_identifiers.csv')
# identifiers = identifiers %>% 
#   rowid_to_column() %>% 
#   rename(order = rowid)
# 
# K4_data = bind_cols(identifiers, 
#                     K4_Qval)
# 
# K4_melted_data = melt(K4_data, 
#                       id.vars = c('order', 
#                                   'population', 
#                                   'individual_id')) %>% 
#   as_tibble() 
# 
# admixture_K4_plot = ggplot(data = K4_melted_data, 
#        aes(x = reorder(order, 
#                        individual_id),
#            y = value, 
#            fill = variable))+
#   geom_bar(stat = "identity", 
#            width = 1)+
#   scale_fill_manual(values = K4_cols)+
#   # scale_fill_manual(values = magma(n = 4))+
#   labs(x = 'Individuals', 
#        y = 'Ancestry proportion')+
#   theme(axis.text.y = element_text(color = 'black'),
#         axis.text.x = element_blank(),
#         axis.title.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         ## can add xaxis labels if needed
#         # axis.text.x = element_text(angle = 90,
#         #                            hjust = 1,
#         #                            vjust = -0.09,
#         #                            size = 8,
#         #                            color = 'black'),
#         legend.position = 'none')+
#   scale_x_discrete(guide = guide_axis(n.dodge = 5))+
#   scale_y_continuous(expand = c(0,0))
# 
# ggsave('admixture_k4_GTS_MYV_SkR_other.tiff',
#        plot = admixture_K4_plot, 
#        dpi = 'retina', 
#        units = 'cm', 
#        width = 35, 
#        height = 10)
# 

##
# pcadapt analysis --------------------------------------------------------

stickle_pcadapt = read.pcadapt('stickleback_maf0.05_ldpruned_filtered.bed', 
                               type = 'bed')

stickle_pca = pcadapt::pcadapt(stickle_pcadapt, 
                               K = 10, 
                               method = 'mahalanobis', 
                               min.maf = 0.01)

plot(stickle_pca, 
     option = 'screeplot', 
     K = 10)

## looks like k = 4 is the best k based on pcadapt

# plot(stickle_pca, option = "scores")

summary(stickle_pca)
## percent variation explained by the three significant axes
stickle_pca$singular.values
sum(stickle_pca$singular.values)
(0.2810947/1.463374)*100
(0.2188855/1.463374)*100
(0.1589348/1.463374)*100
## singular values should have already been square root transformed

# (sqrt(0.2810947)/1.463374)*100
# (sqrt(0.2188855)/1.463374)*100
# (sqrt(0.1589348)/1.463374)*100


stickle_pca$scores

stickle_scores = as_tibble(stickle_pca$scores) %>%
  rename(PC1 = 1,
         PC2 = 2,
         PC3 = 3,
         PC4 = 4,
         PC5 = 5) %>%
  dplyr::select(PC1,
         PC2,
         PC3,
         PC4,
         PC5) %>%
  write_csv('pcadapt_stickle_pc_scores.csv')

stickle_plot = read_csv('pcadapt_stickle_pc_scores.csv')

identifiers = read_csv('stickleback_identifiers.csv')
# identifiers = read_table2('stickleback_maf0.05-ldpruned_nomissing.fam', 
#          col_names = F) %>% 
#   dplyr::select(X1, 
#                 X2) %>% 
#   rename(population = X1, 
#          individual_id = X2) %>% 
#   separate(col = population, 
#            into = c('garbage', 
#                 'population'), 
#            sep = '-') %>% 
#   separate(col = garbage, 
#             into = c('garbage', 
#                    'temp_pop'), 
#             sep = 'Sample_')
# ## going to need to separate into two data frames and 
# ## then bring them back together due to sample name issues
# 
# id_one = identifiers %>% 
#   slice(1:86) %>% 
#   select(population, 
#          individual_id) %>% 
#   separate(col = population, 
#            into = c('population', 
#                     'garbage'), 
#            sep = '_') %>% 
#   slice(1:16) %>% 
#   # separate(individual_id, 
#   #          into = c('garbage', 
#   #                   'individual_id'), 
#   #          sep = 'Sample_') %>% View() 
#   select(population, 
#          individual_id)
# 
# id_three = identifiers %>% 
#   slice(17:86) %>% 
#   select(population, 
#          individual_id) %>% 
#   separate(col = population, 
#            into = c('population', 
#                     'garbage'), 
#            sep = '_') %>% 
#   # slice(1:16) %>% 
#   separate(individual_id,
#            into = c('garbage',
#                     'individual_id'),
#            sep = 'Sample_') %>%
#   select(population, 
#          individual_id)
# 
# id_two = identifiers %>% 
#   slice(87:109) %>% 
#   select(temp_pop, 
#          individual_id) %>% 
#   separate(col = temp_pop, 
#            into = c('population', 
#                     'garbage1', 
#                     'garbage2'), 
#            sep = '_') %>% 
#   separate(col = individual_id, 
#            into = c('garbage', 
#                     'individual_id'), 
#            sep = 'Sample_') %>% 
#   select(population, 
#          individual_id)
# 
# identifiers = bind_rows(id_one,
#                         id_three,
#                         id_two)

# identifiers %>% write_csv('stickleback_identifiers.csv')

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


stickle_plot %>% 
  group_by(Location) %>% 
  summarize(n = n())


##
# pca plot ----------------------------------------------------------------

theme_set(theme_bw())

# cold_warm_cols = c('#264653',
#          '#e63946', 
#          '#386641', 
#          '#6a994e', 
#          '#457b9d', 
#          '#e76f51',
#          '#6d6875', 
#          '#c1121f')


location_cols = c('#00798c',
                  '#003d5b',
                  '#edae49',
                  '#d1495b',
                  '#30638e')


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
  labs(x = 'Principal component 1 (19.21%)',
       y = 'Principal component 2 (14.96%)', 
       col = 'Populations')

stickleback_pca


ggsave(file = 'stickleback_pca_cold_warm_11.08.2023.tiff', 
       path = '~/Parsons_Postdoc/Stickleback_Genomic/Figures/', 
       plot = stickleback_pca, 
       dpi = 'retina', 
       units = 'cm', 
       width = 20.0, 
       height = 13)

# Combine PCA and Admixture -----------------------------------------------

pop_structure_plot = stickleback_pca + admixture_k3_plot

ggsave(file = 'stickleback_population_structure.tiff', 
       path = '~/Parsons_Postdoc/Stickleback_Genomic/Figures/', 
       plot = pop_structure_plot, 
       dpi = 'retina', 
       units = 'cm', 
       width = 50.0, 
       height = 20)

##
# pcadapt outliers --------------------------------------------------------

plot(stickle_pca, 
     option = 'manhattan')

snp_pvalues = stickle_pca$pvalues %>% 
  as_tibble() %>% 
  rename(pvalues = value)

qvalues = qvalue(snp_pvalues$pvalues)
qvalues = qvalues$qvalues 
qvalues = qvalues %>% 
  as_tibble()

snp_vals = bind_cols(snp_pvalues, 
                    qvalues) %>% 
  rename(qvalues = value)

map = read_tsv('stickleback_maf0.05_ldpruned_filtered.map', 
               col_names = F) %>% 
  rename(chromosome = X1, 
         snp = X2, 
         genetic_pos = X3, 
         physical_pos = X4)

final_df = bind_cols(snp_vals, 
          map)

final_df %>% write_csv('pcadapt_qvalues_all_snps.csv')

## calculate cumulative base pair per chromosome
dist_cal = final_df %>% 
  group_by(chromosome) %>% 
  summarise(chr_len = max(physical_pos)) %>% 
  mutate(total = cumsum(chr_len)-chr_len) %>% 
  dplyr::select(-chr_len) %>% 
  left_join(final_df, ., by = c('chromosome'='chromosome')) %>%
  arrange(chromosome, 
          physical_pos) %>% 
  mutate(BPcum = physical_pos + total) 

## calculate the center of the chromosome
axisdf = dist_cal %>% 
  group_by(chromosome) %>% 
  summarize(center=(max(BPcum) + min(BPcum))/2 )  


## Get the neutral snps
# non_outs = dist_cal %>% 
#   filter(qvalues >= 0.05) %>% 
#   mutate
# ## Get the outliers
# outs = dist_cal %>% 
#   filter(qvalues < 0.05)

non_outs = dist_cal %>% 
  filter(qvalues >= 0.05) %>% 
  mutate(logqval = -log(qvalues))
## Get the outliers
outs = dist_cal %>% 
  filter(qvalues < 0.05) %>% 
  mutate(logqval = -log(qvalues))

# outs %>% write_csv('pcadapt_outliers_q0.05.csv')


pcadapt_man = ggplot(non_outs, 
                       aes(x = BPcum, 
                           y = logqval))+
  # plot the non outliers in grey
  geom_point(aes(color = as.factor(chromosome)), 
             alpha = 0.8, 
             size = 1.3)+
  ## alternate colors per chromosome
  scale_color_manual(values = rep(c("grey", "dimgrey"), 39))+
  ## plot the outliers on top of everything
  ## currently digging this hot pink colour
  geom_point(data = outs,
             col = '#8C0F26',
             alpha=0.8, 
             size=1.3)+
  scale_x_continuous(label = axisdf$chromosome, 
                     breaks = axisdf$center)+
  scale_y_continuous(expand = c(0, 0))+
  # scale_y_reverse(expand = c(0, 0))+
  # remove space between plot area and x axis
  labs(x = 'Cumulative base pair', 
       y = '-log10(q-values)')+
  theme(legend.position="none",
        # panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(), 
        axis.text.x = element_text(size = 9, 
                                  angle = 90), 
        axis.title = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12))

pcadapt_man
## ggsave that plot

ggsave(file = 'stickleback_manhattan_plot.tiff', 
       path = 'C:/Stickleback_Genomic/Figures/', 
       plot = pcadapt_man, 
       dpi = 'retina', 
       units = 'cm', 
       width = 20, 
       height = 20)





# pcadapt MYV -------------------------------------------------------------

MYVN_pcadapt = read.pcadapt('MYVN_SNPS.bed', 
                                 type = 'bed')

MYVN_pca = pcadapt::pcadapt(MYVN_pcadapt, 
                               K = 10, 
                               method = 'mahalanobis', 
                               min.maf = 0.01)

plot(MYVN_pca, 
     option = 'screeplot', 
     K = 10)


# plot(stickle_pca, option = "scores")

summary(MYVN_pca)
## percent variation explained by the three significant axes
MYVN_pca$singular.values
sum(MYVN_pca$singular.values)
(0.3486600/2.250287)*100
(0.2720142/2.250287)*100
(0.2398655/2.250287)*100
## singular values should have already been square root transformed

# (sqrt(0.2810947)/1.463374)*100
# (sqrt(0.2188855)/1.463374)*100
# (sqrt(0.1589348)/1.463374)*100


MYVN_scores = as_tibble(MYVN_pca$scores) %>%
  rename(PC1 = 1,
         PC2 = 2,
         PC3 = 3,
         PC4 = 4,
         PC5 = 5) %>%
  dplyr::select(PC1,
                PC2,
                PC3,
                PC4,
                PC5) %>%
  write_csv('MYVN_pc_scores.csv')

MYVN_plot = read_csv('MYVN_pc_scores.csv')

identifiers = read_csv('stickleback_identifiers.csv') %>% 
  filter(population %in% c('MYVC', 
                           'MYVW'))


MYVN_plot = bind_cols(identifiers, 
                         MYVN_plot)

##

# pcadapt outliers MYV --------------------------------------------------------

plot(MYVN_pca, 
     option = 'manhattan')

snp_pvalues = MYVN_pca$pvalues %>% 
  as_tibble() %>% 
  rename(pvalues = value)

qvalues = qvalue(snp_pvalues$pvalues)
qvalues = qvalues$qvalues 
qvalues = qvalues %>% 
  as_tibble()

snp_vals = bind_cols(snp_pvalues, 
                     qvalues) %>% 
  rename(qvalues = value)

map = read_tsv('stickleback_maf0.05_ldpruned_filtered.map', 
               col_names = F) %>% 
  rename(chromosome = X1, 
         snp = X2, 
         genetic_pos = X3, 
         physical_pos = X4)

final_df = bind_cols(snp_vals, 
                     map)


## calculate cumulative base pair per chromosome
dist_cal = final_df %>% 
  group_by(chromosome) %>% 
  summarise(chr_len = max(physical_pos)) %>% 
  mutate(total = cumsum(chr_len)-chr_len) %>% 
  dplyr::select(-chr_len) %>% 
  left_join(final_df, ., by = c('chromosome'='chromosome')) %>%
  arrange(chromosome, 
          physical_pos) %>% 
  mutate(BPcum = physical_pos + total) 

## calculate the center of the chromosome
axisdf = dist_cal %>% 
  group_by(chromosome) %>% 
  summarize(center=(max(BPcum) + min(BPcum))/2 )  


## Get the neutral snps
# non_outs = dist_cal %>% 
#   filter(qvalues >= 0.05) %>% 
#   mutate
# ## Get the outliers
# outs = dist_cal %>% 
#   filter(qvalues < 0.05)

non_outs = dist_cal %>% 
  filter(qvalues >= 0.01) %>% 
  mutate(logqval = -log(qvalues))
## Get the outliers
outs = dist_cal %>% 
  filter(qvalues < 0.01) %>% 
  mutate(logqval = -log(qvalues))

# outs %>% write_csv('pcadapt_outliers_q0.05.csv')

location_cols = c('#06d6a0',
                           '#264653',
                           '#219ebc',
                           '#d62828',
                           '#5f0f40')
                           

MYV_pcadapt_man = ggplot(non_outs, 
                     aes(x = BPcum, 
                         y = logqval))+
  # plot the non outliers in grey
  geom_point(aes(color = as.factor(chromosome)), 
             alpha = 0.8, 
             size = 1.3)+
  ## alternate colors per chromosome
  scale_color_manual(values = rep(c("grey", "dimgrey"), 39))+
  ## plot the outliers on top of everything
  ## currently digging this hot pink colour
  geom_point(data = outs,
             col = '#d62828',
             alpha=0.8, 
             size=1.3)+
  scale_x_continuous(label = axisdf$chromosome, 
                     breaks = axisdf$center)+
  scale_y_continuous(expand = c(0, 0))+
  # scale_y_reverse(expand = c(0, 0))+
  # remove space between plot area and x axis
  labs(x = 'Cumulative base pair', 
       y = '-log10(q-values)')+
  theme(legend.position="none",
        # panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(), 
        axis.text.x = element_text(size = 9, 
                                   angle = 90), 
        axis.title = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12))

MYV_pcadapt_man
## ggsave that plot

ggsave(file = 'stickleback_manhattan_plot.tiff', 
       path = 'C:/Stickleback_Genomic/Figures/', 
       plot = pcadapt_man, 
       dpi = 'retina', 
       units = 'cm', 
       width = 20, 
       height = 20)




# snmf analysis -----------------------------------------------------------
setwd('~/Parsons_Postdoc/Stickleback_Genomic/vcf_filter/')

convert = ped2geno('stickleback_maf0.05_ldpruned_filtered.ped',
                   'stickleback_data.geno')


snmf('stickleback_data.geno',
     K = 1:10,
     entropy = T,
     repetitions = 5,
     project = 'new')

project = load.snmfProject("stickleback_data.snmfProject")
# project = load.snmfProject('Charr_Poly.snmfProject')
# 
summary(project)

ce_plot = plot(project,
     cex = 1.2,
     col = "black",
     pch = 19)

# K=23 has the lowst cross entropy coefficient
ce = cross.entropy(project, K = 5)
best_ce = which.min(ce)

qmatrix = Q(project, 
            K = 5, 
            run = best_ce)

qmatrix %>% 
  as_tibble()
# head(qmatrix)
# dim(qmatrix)
## This will show the clusters each individual was assigned to
## The order is the same as the .ped file!!
## Take the first two columns of the ped file to group by population
apply(qmatrix, 1, which.max) %>%
  as_tibble() %>%
  dplyr::rename(Genetic_group = value) %>%
  write_tsv('stickleback_snmf_k5.txt')

## See if we can make this plot look nice 
## Needs colour and a legend
## internet is out due to a dickhead electrician
## shitty base R plot

location_cols = c('#5f0f40',
                  '#d62828',
                  '#06d6a0',
                  '#264653',
                  '#219ebc')
plot = barchart(project,
               K = 5, 
               run = best_ce, 
               border = NA, 
               space = 0, 
               col = location_cols,
               # sort.by.Q = T,
               # xlab = "Individuals",
               ylab = "Ancestry proportions")
axis(1, at = 1:length(plot$order), 
     labels = plot$order, 
     las = 3, 
     cex.axis = .4)



q_values_k5 = as_tibble(qmatrix) %>%
  dplyr::rename(Q1 = V1,
                Q2 = V2,
                Q3 = V3,
                Q4 = V4, 
                Q5 = V5) 

plot_order = plot$order %>% 
  as_tibble()
  
bind_cols(plot_order, 
          q_values_k5) %>% 
  rename(order = 1) %>% 
  write_csv('stickleback_snmf_qvalues_k5.csv')


snmf_data = read_csv('stickleback_snmf_qvalues_k5.csv')

identifiers %>% 
  group_by(population) %>% 
  summarize(n = n())

identifiers = read_csv('stickleback_identifiers.csv')
identifiers = identifiers %>% 
  rowid_to_column() %>% 
  rename(order = rowid)

View(identifiers)
# View(identifiers)
# View(plot_order)
snmf_data = left_join(snmf_data, 
          identifiers, 
          by = 'order') 
# View(snmf_data)
snmf_melted = melt(snmf_data, 
                   id.vars = c('order',
                               'population', 
                               'individual_id')) %>% 
  as_tibble() %>% 
  group_by(population) %>% 
  arrange(population)

snmf_melted$order = as.character(snmf_melted$order)

## need colour scheme
## real
location_cols = c('#264653',
                  '#5f0f40',
                  '#06d6a0',
                  '#219ebc',
                  '#d62828')

# snmf_melted %>%
#   select(variable) %>%
#   distinct()
data_group1 = snmf_melted %>% 
  ungroup() %>% 
  arrange(individual_id) %>% 
  # View()
  slice(1:430) %>% 
  # View()
# write_csv('snmf_stickleback_k5_data.csv')
  separate('individual_id', 
           c('garbage',
             'individual_id'), 
           sep = '-') %>% 
  select(-garbage)

data_group2 = snmf_melted %>% 
  ungroup() %>% 
  arrange(individual_id) %>% 
  slice(430:545) 

snmf_melted_fixed = bind_rows(data_group1, 
          data_group2) %>% 
  arrange(order)

View(snmf_melted_fixed)

## snmf plot
snmf_plot = ggplot(data = snmf_melted_fixed, 
                   aes(x = fct_inorder(individual_id),
                       y = value, 
                       fill = variable, 
                       group = population))+
  geom_bar(stat = "identity", 
           width = 1)+
  scale_fill_manual(values = location_cols)+
  # scale_fill_manual(values = magma(n = 4))+
  labs(x = 'Individuals', 
       y = 'Ancestry proportion')+
  theme(axis.text.y = element_text(color = 'black'),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        ## can add xaxis labels if needed
        # axis.text.x = element_text(angle = 90,
        #                            hjust = 1,
        #                            vjust = -0.09,
        #                            size = 10,
        #                            color = 'black'),
        legend.position = 'none')+
  scale_x_discrete(guide = guide_axis(n.dodge = 5))+
  scale_y_continuous(expand = c(0,0))

snmf_plot

ggsave(file = 'stickleback_snmf.tiff', 
       path = 'C:/Stickleback_Genomic/Figures/', 
       plot = snmf_plot, 
       dpi = 'retina', 
       units = 'cm', 
       width = 30, 
       height = 8)




# PCA Umap --------------------------------------------------------------------

stickle_pca = read.pcadapt('stickleback_maf0.05_ldpruned_filtered.bed', 
                               type = 'bed')

## First we need to compute the first 50 pc axes
## based on the umap paper
stickle_pca = pcadapt::pcadapt(stickle_pca, 
                               K = 50, 
                               method = 'mahalanobis', 
                               min.maf = 0.01)

plot(stickle_pca, 
     option = 'screeplot', 
     K = 50)

pcscores = as_tibble(stickle_pca$scores) 

pc_data = bind_cols(identifiers, 
                    pcscores)

# pc_data %>% write_csv('stickleback_pcscores_umap.csv')

## now we can run umap on the pc scores
set.seed(666)

umap_fit = pc_data %>%
  dplyr::select(V1:V50) %>% 
  scale() %>% 
  umap()


##
# Fst set up Plink ---------------------------------------------------------------

ped_test = read_table2('stickleback_maf0.05_ldpruned_filtered.ped', 
            col_names = F) 

# read_tsv('stickleback_maf0.05_ldpruned_filtered.map', 
#                col_names = F) %>% 
#   filter(X1 == 'chr_XXI') %>% 
#   filter(X4 >= 9963830, 
#          X4 <= 11574024) %>% 
#   write_tsv('stickleback_afvaper_500snps.map')
  

identifiers = read_csv('stickleback_identifiers.csv')
identifiers = mutate(.data = identifiers,
       type = as.factor(case_when(
         population == 'ASHNC' ~ 'Cold',
         population == 'ASHNW' ~ 'Warm',
         population == 'CSWY' ~ 'Manmade',
         population == 'GTS' ~ 'Thermal',
         population == 'MYVC' ~ 'Cold',
         population == 'MYVW' ~ 'Warm',
         population == 'SKRC' ~ 'Cold',
         population == 'SKRW' ~ 'Warm')))

ped_ids = read_table2('stickleback_maf0.05_ldpruned_filtered.fam', 
                      col_names = F) %>%
  dplyr::select(X1,
                X2)
ped_ids = bind_cols(ped_ids, 
                    identifiers)

ped_ids %>% 
  select(individual_id) %>% 
  View()

## Need to split based on each comparison 
## ASHW vs ASHC
## MYVW vs MYVC
## SKRW vs SKRC
## GTS vs CSWY

# ped_ids %>% 
#   filter(type %in% c('Warm', 
#                          'Cold')) %>%
#   select(X1, 
#          X2, 
#          type) %>% 
#   rename(`#population` = 1, 
#          individual_id = 2) %>% 
#   write_tsv('Warm_cold_Fst_grouping.txt')

## Holy fuck!! Make sure to use the actual family and individual
## identifiers in the fucking ped file. WOW

 

## Need to make a ped and map file for each of these comparisons
## the ped file is waaaay to big to open in R
## use the --keep or --keep-fam flags in plink to filter the 
## populations out. 
## the --keep file needs to be a text file with family and individual
## identifiers

# temp = read_tsv('temp.env', 
#                 col_names = 'temp')
# 
# ped_ids = bind_cols(ped_ids, 
#                     temp)

# ped_ids %>% 
#   filter(population %in% c('GTS', 
#                            'CSWY')) %>% 
#   dplyr::select(temp) %>% 
#   write_tsv('GTS_CSWY_temp_var.env', 
#             col_names = F)

ped_ids %>% 
  filter(population %in% c('MYVC', 
                           'MYVW', 
                           'ASHNC', 
                           'ASHNW', 
                           'SKRC', 
                           'SKRW')) %>% 
  # filter(population %in% c('GTS', 
  #                          'CSWY')) %>% 
  # filter(type %in% c('Warm', 
  #                          'Cold')) %>%
  dplyr::select(X1, 
                X2) %>% 
  # rename(`#population` = population, 
         # individual_ID = X1) %>% 
  write_tsv('Warm_Cold_No_GTS_CSWY_keep.txt', 
            col_names = F)

##
# FST analysis ------------------------------------------------------------

ASHN_Fst = read_tsv('ASHN_Fst_values.fst') %>% 
  na.omit() %>%  ##pull out na's
  mutate(FST_zero = if_else(FST < 0, 0, FST))
MYV_Fst = read_tsv('MYV_Fst_values.fst')%>% 
  na.omit() %>%  ##pull out na's
  mutate(FST_zero = if_else(FST < 0, 0, FST))
SKR_Fst = read_tsv('SKR_Fst_values.fst')%>% 
  na.omit() %>%  ##pull out na's
  mutate(FST_zero = if_else(FST < 0, 0, FST))
GTS_CSWY_Fst = read_tsv('GTS_CSWY_Fst_values.fst')%>% 
  na.omit() %>%  ##pull out na's
  mutate(FST_zero = if_else(FST < 0, 0, FST))
WC_Fst = read_tsv('Warm_Cold_Fst.fst') %>% 
  na.omit() %>% 
  mutate(FST_zero = if_else(FST < 0, 0, FST))

##
# FST outliers ------------------------------------------------------------

## Top 5% of the distribution of Fst values
# WC_top_dist = WC_Fst[WC_Fst$FST_zero > quantile(WC_Fst$FST_zero, 
#                                                       prob = 1-5/100),]
## The top 5% might not be restrictive enough to be honest
## What does the top 0.5% look like? 
WC_top_dist = WC_Fst[WC_Fst$FST_zero > quantile(WC_Fst$FST_zero, 
                                                prob = 1-0.5/100),]
# ASHN_top_dist %>% write_csv('ASHN_FST_Outliers.csv')
WC_top_dist %>% 
  # group_by(CHR) %>% 
  summarise(min_fst = min(FST_zero), 
            max_fst = max(FST_zero), 
            mean_fst = mean(FST_zero)) 

WC_Top_Dawg_Outs = WC_top_dist %>% 
  filter(FST_zero >= 0.0982) %>% 
  arrange(CHR, 
          POS) 

# WC_Top_Dawg_Outs %>%
#   dplyr::select(CHR,
#          SNP,
#          POS,
#          FST_zero) %>%
#   dplyr::rename(FST = FST_zero) %>%
#   write_csv('WC_FST_Outliers_WARM_COLD_COMBO.csv')

## snps that are the top 5% fst distribution
# ASHN_top_dist = ASHN_Fst[ASHN_Fst$FST_zero > quantile(ASHN_Fst$FST_zero, 
#                                     prob = 1-5/100),]

ASHN_top_dist = ASHN_Fst[ASHN_Fst$FST_zero > quantile(ASHN_Fst$FST_zero, 
                                                      prob = 1-0.5/100),]

# ASHN_top_dist %>% write_csv('ASHN_FST_Outliers.csv')
ASHN_top_dist %>% 
  # group_by(CHR) %>% 
  summarise(min_fst = min(FST_zero), 
            max_fst = max(FST_zero), 
            mean_fst = mean(FST_zero)) 

ASHN_Top_Dawg_Outs = ASHN_top_dist %>% 
  filter(FST_zero >= 0.235) %>% 
  arrange(CHR, 
          POS) 
ASHN_Top_Dawg_Outs %>%
  dplyr::select(CHR,
                SNP,
                POS,
                FST_zero) %>%
  dplyr::rename(FST = FST_zero) %>%
  write_csv('FST_Outliers_ASHN.csv')

# MYV_top_dist = MYV_Fst[MYV_Fst$FST_zero > quantile(MYV_Fst$FST_zero, 
#                                                       prob = 1-5/100),]

MYV_top_dist = MYV_Fst[MYV_Fst$FST_zero > quantile(MYV_Fst$FST_zero, 
                                                   prob = 1-0.5/100),]

# MYV_top_dist %>% write_csv('MYV_FST_Outliers.csv')
MYV_top_dist %>% 
  # group_by(CHR) %>% 
  summarise(min_fst = min(FST_zero), 
            max_fst = max(FST_zero), 
            mean_fst = mean(FST_zero))

MYV_Top_Dawg_Outs = MYV_top_dist %>% 
  filter(FST_zero >= 0.238) %>% 
  arrange(CHR, 
          POS)

# MYV_Top_Dawg_Outs %>% 
#   dplyr::select(CHR, 
#                 SNP, 
#                 POS, 
#                 FST_zero) %>% 
#   dplyr::rename(FST = FST_zero) %>% 
#   write_csv('FST_Outliers_MYV.csv')
# SKR_top_dist = SKR_Fst[SKR_Fst$FST_zero > quantile(SKR_Fst$FST_zero, 
#                                                    prob = 1-5/100),]
SKR_top_dist = SKR_Fst[SKR_Fst$FST_zero > quantile(SKR_Fst$FST_zero, 
                                                   prob = 1-0.5/100),]
# SKR_top_dist %>% write_csv('SKR_FST_Outliers.csv')

SKR_top_dist %>% 
  # group_by(CHR) %>% 
  summarise(min_fst = min(FST_zero), 
            max_fst = max(FST_zero), 
            mean_fst = mean(FST_zero))

SKR_Top_Dawg_Outs = SKR_top_dist %>% 
  filter(FST_zero >= 0.306) %>% 
  arrange(CHR, 
          POS)
# SKR_Top_Dawg_Outs %>% 
#   dplyr::select(CHR, 
#                 SNP, 
#                 POS, 
#                 FST_zero) %>% 
#   dplyr::rename(FST = FST_zero) %>% 
#   write_csv('FST_Outliers_SKR.csv')

# GTS_CSWY_top_dist = GTS_CSWY_Fst[GTS_CSWY_Fst$FST_zero > quantile(GTS_CSWY_Fst$FST_zero, 
#                                                    prob = 1-5/100),]

GTS_CSWY_top_dist = GTS_CSWY_Fst[GTS_CSWY_Fst$FST_zero > quantile(GTS_CSWY_Fst$FST_zero, 
                                                                  prob = 1-0.5/100),]

# GTS_CSWY_top_dist %>% write_csv('GTS_CSWY_FST_Outliers.csv')

GTS_CSWY_top_dist %>% 
  # group_by(CHR) %>% 
  summarise(min_fst = min(FST_zero), 
            max_fst = max(FST_zero), 
            mean_fst = mean(FST_zero))

GTS_CSWY_Top_Dawg_Outs = GTS_CSWY_top_dist %>% 
  filter(FST_zero >= 0.706) %>% 
  arrange(CHR, 
          POS)
# GTS_CSWY_Top_Dawg_Outs %>% 
#   dplyr::select(CHR, 
#                 SNP, 
#                 POS, 
#                 FST_zero) %>% 
#   dplyr::rename(FST = FST_zero) %>% 
#   write_csv('FST_Outliers_GTS_CSWY.csv')
##### TOP DAWG OUTS Fst manhattan set up ----------------------------------------------------

Fst_manhatan_format(ASHN_Fst,
                    ASHN_Top_Dawg_Outs) %>%
  write_csv('ASHN_TOP_DAWG_Fst_clean.csv')

Fst_manhatan_format(MYV_Fst,
                    MYV_Top_Dawg_Outs) %>% 
  write_csv('MYV_TOP_DAWG_Fst_clean.csv')

Fst_manhatan_format(SKR_Fst,
                    SKR_Top_Dawg_Outs) %>% 
  write_csv('SKR_TOP_DAWG_Fst_clean.csv')

Fst_manhatan_format(GTS_CSWY_Fst,
                    GTS_CSWY_Top_Dawg_Outs) %>% 
  write_csv('GTS_CSWY_TOP_DAWG_Fst_clean.csv')

# WC_Top_Dawg_Outs

Fst_manhatan_format(WC_Fst,
                    WC_Top_Dawg_Outs) %>% 
  write_csv('WC_TOP_DAWG_Fst_clean.csv')
# TOP DAWG OUTS FST outlier manhattan plot ----------------------------------------------

ASHN_Fst_clean = read_csv('ASHN_TOP_DAWG_Fst_clean.csv') %>% 
  stickle_CHR_reorder() %>% 
  dist_cal()
MYV_Fst_clean = read_csv('MYV_TOP_DAWG_Fst_clean.csv') %>% 
  stickle_CHR_reorder() %>% 
  dist_cal()
SKR_Fst_clean = read_csv('SKR_TOP_DAWG_Fst_clean.csv') %>% 
  stickle_CHR_reorder() %>% 
  dist_cal()
GTS_CSWY_Fst_clean = read_csv('GTS_CSWY_TOP_DAWG_Fst_clean.csv') %>% 
  stickle_CHR_reorder() %>% 
  dist_cal()

## calculate the center of the chromosome
axisdf_ASHN = axis_df(ASHN_Fst_clean)
axisdf_MYV = axis_df(MYV_Fst_clean)
axisdf_SKR = axis_df(SKR_Fst_clean)
axisdf_GTS = axis_df(GTS_CSWY_Fst_clean)

non_outs = 
  ASHN_Fst_clean %>%
  # MYV_Fst_clean %>%
  # SKR_Fst_clean %>%
  # GTS_CSWY_Fst_clean %>%
  filter(value == 'Neutral') 

## Get the outliers
outs = 
  ASHN_Fst_clean %>%
  # MYV_Fst_clean %>%
  # SKR_Fst_clean %>%
  # GTS_CSWY_Fst_clean %>%
  filter(value == 'Outlier') 

## ASHN colour = #00798c
## MYV colour = #d1495b
## SKR colour = #30638e
## GTS_CSWY colour = #edae49
## WC colour = #439a86


ASHN_Fst_manhattan = Fst_manhattan(non_outs = non_outs, 
                                   outs = outs, 
                                   axisdf = axisdf_ASHN, 
                                   xval = BPcum, 
                                   yval = FST_zero, 
                                   chr = non_outs$CHR, 
                                   out_col = '#00798c', 
                                   plot_letter = 'A) Áshildarholtsvatn')

MYV_Fst_manhattan = Fst_manhattan(non_outs = non_outs, 
                                  outs = outs, 
                                  axisdf = axisdf_MYV, 
                                  xval = BPcum, 
                                  yval = FST_zero, 
                                  chr = non_outs$CHR, 
                                  out_col = '#d1495b', 
                                  plot_letter = 'B) Mývatn')

SKR_Fst_manhattan = Fst_manhattan(non_outs = non_outs, 
                                  outs = outs, 
                                  axisdf = axisdf_SKR, 
                                  xval = BPcum, 
                                  yval = FST_zero, 
                                  chr = non_outs$CHR, 
                                  out_col = '#30638e', 
                                  plot_letter = 'C) Sauðárkrókur')

GTS_CSWY_Fst_manhattan = Fst_manhattan(non_outs = non_outs,
                                       outs = outs, 
                                       axisdf = axisdf_GTS, 
                                       xval = BPcum, 
                                       yval = FST_zero, 
                                       chr = non_outs$CHR, 
                                       out_col = '#edae49', 
                                       plot_letter = 'D) Grettislaug-Garðsvatn')


TOP_DAWG_PLOT = (ASHN_Fst_manhattan|MYV_Fst_manhattan)/(SKR_Fst_manhattan|GTS_CSWY_Fst_manhattan)


## ggsave that plot

ggsave(file = 'stickleback_TOP_DAWG_0.5%_FST_manhattan_plot.tiff', 
       path = '~/Parsons_Postdoc/Stickleback_Genomic/Figures/', 
       plot = TOP_DAWG_PLOT, 
       dpi = 'retina', 
       units = 'cm', 
       width = 40, 
       height = 20)


# Fst violin plots --------------------------------------------------------

library(qqman)

manhattan(gwasResults, chr="CHR", bp="BP", snp="SNP", p="P" )

ASHN_Fst_clean

test = Chr_convert2(ASHN_Fst_clean)

test$chr_num = as.numeric(test$chr_num)

manhattan(test, 
          chr = 'chr_num', 
          bp = 'POS', 
          snp = 'SNP', 
          p = 'FST_zero')

MYV_Fst_clean %>% 
  mutate(test_name = 'Group') %>% 
  ggplot()+
  geom_dotplot(aes(x = test_name, 
                   y = FST_zero))
  geom_violin(aes(x = test_name, 
                  y = FST_zero))

MYV_Fst_clean %>% 
  filter(FST_zero < 0.002)

MYV_Fst_clean %>% 
  filter(FST_zero > 0.4, 
         FST_zero < 0.5)

ASHN_Fst_clean %>% 
  filter(FST_zero > 0.0022, 
         FST_zero < 0.1)
ASHN_Fst_clean %>% 
  filter(FST_zero > 0.1, 
         FST_zero < 0.2)
ASHN_Fst_clean %>% 
  filter(FST_zero > 0.2, 
         FST_zero < 0.3)

ASHN_Fst_clean %>% 
  filter(FST_zero > 0.3, 
         FST_zero < 0.4)
ASHN_Fst_clean %>% 
  filter(FST_zero > 0.4, 
         FST_zero < 0.5)
ASHN_Fst_clean %>% 
  filter(FST_zero > 0.5)

  
non_outs = 
  # ASHN_Fst_clean %>%
  # MYV_Fst_clean %>%
  # SKR_Fst_clean %>%
  GTS_CSWY_Fst_clean %>%
  filter(value == 'Neutral') 

## Get the outliers
outs = 
  # ASHN_Fst_clean %>%
  # MYV_Fst_clean %>%
  # SKR_Fst_clean %>%
  GTS_CSWY_Fst_clean %>%
  filter(value == 'Outlier') 


##
# FST distribution plots --------------------------------------------------

location_cols = c('#00798c',
                  '#003d5b',
                  '#edae49',
                  '#d1495b',
                  '#30638e')

ASHN_Fst_dist_plot = ASHN_Fst %>% 
  ggplot()+
  geom_density(aes(x = FST_zero), 
               col = '#00798c', 
               fill = '#00798c')+
  geom_density(data = ASHN_top_dist, 
               aes(x = FST_zero),
               col = '#00798c',
               fill = '#00798c')+
  labs(x = 'Fst', 
       y = 'Density', 
       title = 'A)')+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12))

MYV_Fst_dist_plot = MYV_Fst %>% 
  ggplot()+
  geom_density(aes(x = FST_zero), 
               col = '#d1495b', 
               fill = '#d1495b')+
  geom_density(data = MYV_top_dist, 
               aes(x = FST_zero),
               col = '#d1495b',
               fill = '#d1495b')+
  labs(x = 'Fst', 
       y = 'Density', 
       title = 'B)')+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12))


SKR_Fst_dist_plot = SKR_Fst %>% 
  ggplot()+
  geom_density(aes(x = FST_zero), 
               col = '#30638e', 
               fill = '#30638e')+
  geom_density(data = SKR_top_dist, 
               aes(x = FST_zero),
               col = '#30638e',
               fill = '#30638e')+
  labs(x = 'Fst', 
       y = 'Density', 
       title = 'C)')+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12))


GTS_CSWY_Fst_dist_plot = GTS_CSWY_Fst %>% 
  ggplot()+
  geom_density(aes(x = FST_zero), 
               col = '#edae49', 
               fill = '#edae49')+
  geom_density(data = GTS_CSWY_top_dist, 
               aes(x = FST_zero),
               col = '#edae49',
               fill = '#edae49')+
  labs(x = 'Fst', 
       y = 'Density', 
       title = 'D)')+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12))


WC_Fst_dist_plot = WC_Fst %>% 
  ggplot()+
  geom_density(aes(x = FST_zero), 
               col = '#439a86', 
               fill = '#439a86')+
  geom_density(data = WC_top_dist, 
               aes(x = FST_zero),
               col = '#439a86',
               fill = '#439a86')+
  labs(x = 'Fst', 
       y = 'Density', 
       title = 'E)')+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12))


Fst_dist_combo = (ASHN_Fst_dist_plot|MYV_Fst_dist_plot)/(SKR_Fst_dist_plot|GTS_CSWY_Fst_dist_plot) | WC_Fst_dist_plot 

ggsave(file = 'stickleback_0.5%_FST_Distribution_plots.tiff', 
       path = 'C:/Stickleback_Genomic/Figures/', 
       plot = Fst_dist_combo, 
       dpi = 'retina', 
       units = 'cm', 
       width = 20.0, 
       height = 15)

# Fst distribution per chromosome plot ------------------------------------
location_cols = c('#00798c',
                  '#003d5b',
                  '#edae49',
                  '#d1495b',
                  '#30638e')

ASHN_per_chrom = ASHN_Fst %>% 
  ggplot()+
  geom_density(aes(x = FST_zero), 
               col = '#00798c', 
               fill = '#00798c')+
  geom_density(data = ASHN_top_dist, 
               aes(x = FST_zero),
               col = '#00798c',
               fill = '#00798c')+
  facet_grid(~CHR)+
  labs(x = 'Fst', 
       y = 'Density', 
       title = 'A)')+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(),
        strip.background = element_rect(fill = 'white'), 
        strip.text = element_text(color = 'black'))


MYV_per_chrom = MYV_Fst %>% 
  ggplot()+
  geom_density(aes(x = FST_zero), 
               col = '#d1495b', 
               fill = '#d1495b')+
  geom_density(data = MYV_top_dist, 
               aes(x = FST_zero),
               col = '#d1495b',
               fill = '#d1495b')+
  facet_grid(~CHR)+
  labs(x = 'Fst', 
       y = 'Density', 
       title = 'B)')+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(),
        strip.background = element_rect(fill = 'white'), 
        strip.text = element_text(color = 'black'))

SKR_per_chrom = SKR_Fst %>% 
  ggplot()+
  geom_density(aes(x = FST_zero), 
               col = '#30638e', 
               fill = '#30638e')+
  geom_density(data = SKR_top_dist, 
               aes(x = FST_zero),
               col = '#30638e',
               fill = '#30638e')+
  facet_grid(~CHR)+
  labs(x = 'Fst', 
       y = 'Density', 
       title = 'C)')+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(),
        strip.background = element_rect(fill = 'white'), 
        strip.text = element_text(color = 'black'))

GTS_CSWY_per_chrom = GTS_CSWY_Fst %>% 
  ggplot()+
  geom_density(aes(x = FST_zero), 
               col = '#edae49', 
               fill = '#edae49')+
  geom_density(data = GTS_CSWY_top_dist, 
               aes(x = FST_zero),
               col = '#000000',
               fill = '#000000')+
  facet_grid(~CHR)+
  labs(x = 'Fst', 
       y = 'Density', 
       title = 'D)')+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(),
        strip.background = element_rect(fill = 'white'), 
        strip.text = element_text(color = 'black'))

WC_per_chrom = WC_Fst %>% 
  ggplot()+
  geom_density(aes(x = FST_zero), 
               col = '#439a86', 
               fill = '#439a86')+
  geom_density(data = WC_top_dist, 
               aes(x = FST_zero),
               col = '#439a86',
               fill = '#439a86')+
  facet_grid(~CHR)+
  labs(x = 'Fst', 
       y = 'Density', 
       title = 'E)')+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 14), 
        axis.text.x = element_text(size = 8, 
                                   angle = 90), 
        axis.text.y = element_blank(), 
        axis.title.y = element_blank(), 
        axis.ticks = element_blank(),
        strip.background = element_rect(fill = 'white'), 
        strip.text = element_text(color = 'black'))


per_chrom_combo = ASHN_per_chrom/MYV_per_chrom/SKR_per_chrom/GTS_CSWY_per_chrom/WC_per_chrom



ggsave(file = 'stickleback_0.5%_FST_Distribution_per_chrome.tiff', 
       path = 'C:/Stickleback_Genomic/Figures/', 
       plot = per_chrom_combo, 
       dpi = 'retina', 
       units = 'cm', 
       width = 40.0, 
       height = 20)


##
# Fst manhattan set up ----------------------------------------------------


# intersect(ASHN_Fst$SNP, 
#           ASHN_top_dist$SNP) %>% as_tibble()


WC_Fst_clean = Fst_manhatan_format(WC_Fst, 
                                   WC_top_dist) %>% 
  write_csv('WC_Fst_clean.csv')

ASHN_Fst_clean = Fst_manhatan_format(ASHN_Fst, 
                                     ASHN_top_dist) %>%
  write_csv('ASHN_Fst_clean.csv')

MYV_Fst_clean = Fst_manhatan_format(MYV_Fst, 
                                    MYV_top_dist) %>% 
  write_csv('MYV_Fst_clean.csv')

SKR_Fst_clean = Fst_manhatan_format(SKR_Fst, 
                                    SKR_top_dist) %>% 
  write_csv('SKR_Fst_clean.csv')

GTS_CSWY_Fst_clean = Fst_manhatan_format(GTS_CSWY_Fst, 
                                         GTS_CSWY_top_dist) %>% 
  write_csv('GTS_CSWY_Fst_clean.csv')

# FST outlier manhattan plot ----------------------------------------------

ASHN_Fst_clean = read_csv('ASHN_Fst_clean.csv') %>% 
  stickle_CHR_reorder() %>% 
  dist_cal()
MYV_Fst_clean = read_csv('MYV_Fst_clean.csv') %>% 
  stickle_CHR_reorder() %>% 
  dist_cal()
SKR_Fst_clean = read_csv('SKR_Fst_clean.csv') %>% 
  stickle_CHR_reorder() %>% 
  dist_cal()
GTS_CSWY_Fst_clean = read_csv('GTS_CSWY_Fst_clean.csv') %>% 
  stickle_CHR_reorder() %>% 
  dist_cal()
WC_Fst_clean = read_csv('WC_Fst_clean.csv') %>% 
  stickle_CHR_reorder() %>% 
  dist_cal()
  ## calculate the center of the chromosome
axisdf_ASHN = axis_df(ASHN_Fst_clean)
axisdf_MYV = axis_df(MYV_Fst_clean)
axisdf_SKR = axis_df(SKR_Fst_clean)
axisdf_GTS = axis_df(GTS_CSWY_Fst_clean)
axisdf_WC = axis_df(WC_Fst_clean)

non_outs = 
  WC_Fst_clean %>%
  # ASHN_Fst_clean %>%
  # MYV_Fst_clean %>%
  # SKR_Fst_clean %>%
  # GTS_CSWY_Fst_clean %>%
  filter(value == 'Neutral') 

## Get the outliers
outs = 
  WC_Fst_clean %>%
  # ASHN_Fst_clean %>%
  # MYV_Fst_clean %>%
  # SKR_Fst_clean %>%
  # GTS_CSWY_Fst_clean %>%
  filter(value == 'Outlier') 

## ASHN colour = #00798c
## MYV colour = #d1495b
## SKR colour = #30638e
## GTS_CSWY colour = #edae49
## WC colour = #439a86

WC_Fst_manhattan = Fst_manhattan(non_outs = non_outs, 
                                 outs = outs, 
                                 axisdf = axisdf_WC, 
                                 xval = BPcum, 
                                 yval = FST_zero, 
                                 chr = non_outs$CHR, 
                                 out_col = '#439a86', 
                                 plot_letter = 'E) Geothermal-Ambient comparison')


ASHN_Fst_manhattan = Fst_manhattan(non_outs = non_outs, 
                                       outs = outs, 
                                       axisdf = axisdf_ASHN, 
                                       xval = BPcum, 
                                       yval = FST_zero, 
                                       chr = non_outs$CHR, 
                                       out_col = '#00798c', 
                                       plot_letter = 'A) Áshildarholtsvatn')

MYV_Fst_manhattan = Fst_manhattan(non_outs = non_outs, 
                                   outs = outs, 
                                   axisdf = axisdf_MYV, 
                                   xval = BPcum, 
                                   yval = FST_zero, 
                                   chr = non_outs$CHR, 
                                   out_col = '#d1495b', 
                                   plot_letter = 'B) Mývatn')

SKR_Fst_manhattan = Fst_manhattan(non_outs = non_outs, 
                                  outs = outs, 
                                  axisdf = axisdf_SKR, 
                                  xval = BPcum, 
                                  yval = FST_zero, 
                                  chr = non_outs$CHR, 
                                  out_col = '#30638e', 
                                  plot_letter = 'C) Sauðárkrókur')

GTS_CSWY_Fst_manhattan = Fst_manhattan(non_outs = non_outs,
              outs = outs, 
              axisdf = axisdf_GTS, 
              xval = BPcum, 
              yval = FST_zero, 
              chr = non_outs$CHR, 
              out_col = '#edae49', 
              plot_letter = 'D) Grettislaug-Garðsvatn')


Fst_man_combo = (ASHN_Fst_manhattan|MYV_Fst_manhattan)/(SKR_Fst_manhattan|GTS_CSWY_Fst_manhattan)|WC_Fst_manhattan


## ggsave that plot

ggsave(file = 'stickleback_0.5%_FST_manhattan_plot.tiff', 
       path = 'C:/Stickleback_Genomic/Figures/', 
       plot = Fst_man_combo, 
       dpi = 'retina', 
       units = 'cm', 
       width = 40, 
       height = 20)


local_adaptation = (ASHN_Fst_manhattan|MYV_Fst_manhattan)/(SKR_Fst_manhattan|GTS_CSWY_Fst_manhattan)


## ggsave that plot

ggsave(file = 'local_adaptation_stickleback_0.5%_FST_manhattan_plot.tiff', 
       path = '~/Parsons_Postdoc/Stickleback_Genomic/Figures/', 
       plot = local_adaptation, 
       dpi = 'retina', 
       units = 'cm', 
       width = 40, 
       height = 20)



# Outlier overlap ---------------------------------------------------------

## per population FST outliers
ASHN_Fst_clean = read_csv('ASHN_TOP_DAWG_Fst_clean.csv') %>%
  stickle_CHR_reorder() %>%
  dist_cal() %>%
  filter(value == 'Outlier') %>% 
  select(SNP)
MYV_Fst_clean = read_csv('MYV_TOP_DAWG_Fst_clean.csv') %>%
  stickle_CHR_reorder() %>%
  dist_cal()%>%
  filter(value == 'Outlier') %>% 
  select(SNP)
SKR_Fst_clean = read_csv('SKR_TOP_DAWG_Fst_clean.csv') %>%
  stickle_CHR_reorder() %>%
  dist_cal()%>%
  filter(value == 'Outlier') %>% 
  select(SNP)
GTS_CSWY_Fst_clean = read_csv('GTS_CSWY_TOP_DAWG_Fst_clean.csv') %>%
  stickle_CHR_reorder() %>%
  dist_cal()%>%
  filter(value == 'Outlier') %>% 
  select(SNP)

WC_FST_clean = read_csv("WC_FST_Outliers_WARM_COLD_COMBO.csv") %>%
  stickle_CHR_reorder() %>%
  dist_cal()%>%
  filter(value == 'Outlier') %>% 
  select(SNP)

intersect(ASHN_Fst_clean, 
          MYV_Fst_clean)

intersect(ASHN_Fst_clean, 
          SKR_Fst_clean)

intersect(ASHN_Fst_clean, 
          GTS_CSWY_Fst_clean)

intersect(MYV_Fst_clean, 
          SKR_Fst_clean)

intersect(MYV_Fst_clean, 
          GTS_CSWY_Fst_clean)

intersect(SKR_Fst_clean, 
          GTS_CSWY_Fst_clean)

##Common FST outliers
WC_Fst_clean_outs = read_csv('WC_Fst_clean.csv') %>% 
  stickle_CHR_reorder() %>% 
  dist_cal()%>% 
  filter(value == 'Outlier')

# WC_Fst_clean_outs %>% 
#   group_by(CHR) %>% 
#   summarise(n_outs = n()) %>% 
#   View()
# # 
# WC_Fst_clean_all = read_csv('WC_Fst_clean.csv') %>% 
#   stickle_CHR_reorder() %>% 
#   dist_cal()

## pcadapt outliers
common_pcadapt_outliers = read_csv('pcadapt_outliers_q0.05.csv') %>% 
  rename(SNP = snp, 
         CHR = chromosome, 
         POS = physical_pos)

# common_pcadapt_outliers %>% 
#   group_by(CHR) %>% 
#   summarise(n_outlier = n()) %>% 
#   View()


FST_out_pcadapt = inner_join(WC_Fst_clean_outs, 
           common_pcadapt_outliers, 
          by = c('CHR', 
                 'SNP', 
                 'POS')) 
# %>% 
#   group_by(CHR) %>% 
#   summarise(n_outlier = n()) %>% 
#   view()

# FST_all_pcadapt_snps = inner_join(WC_Fst_clean_all, 
#            common_pcadapt_outliers, 
#            by = c('CHR', 
#                   'SNP', 
#                   'POS'))
##LFMM outliers
# lfmm_outliers = read.vcfR("C:/Stickleback_Genomic/vcf_filter/lfmm.SNPs.vcf")

LFMM_Outliers = read_csv('LFMM_Mapped_Full_Data.csv') %>% 
  filter(qvalue <= 0.05) %>% 
  mutate(qvalue_trans = -log10(qvalue))


FST_outs_LFMM = inner_join(WC_Fst_clean_outs, 
           LFMM_Outliers, 
           by = c('CHR', 
                  'SNP', 
                  'POS'))

LFMM_pcadapt = inner_join(common_pcadapt_outliers, 
                           LFMM_Outliers, 
                           by = c('CHR', 
                                  'SNP', 
                                  'POS'))

# View(LFMM_pcadapt)

## Damn only 33 outliers from all three analyses
Three_analysis_outs = inner_join(FST_outs_LFMM, 
           common_pcadapt_outliers,
           by = c('CHR', 
                  'SNP', 
                  'POS'))
## two common outliers on chromosome 8 and 19

Three_analysis_outs %>% 
  select(CHR, 
         SNP, 
         POS, 
         FST_zero) %>% 
  arrange(CHR) %>% 
  filter(FST_zero >= 0.04) %>% 
  # group_by(CHR) %>% 
  # summarise(n = n()) %>% 
  View()

## NEED TO MAKE SURE THE CHROMOSOME ORDER FROM AFVAPER LINES 
## UP TO THESE CHROMOSOMES

# Three_analysis_outs %>%
#   group_by(CHR) %>%
#   summarise(n_outs = n(),
#             Fst_chr = mean(FST_zero)) %>%
#   View()

# Three_analysis_outs %>% write_csv('Outliers_FST_pcadapt_LFMM.csv')

# %>% 
#   group_by(CHR) %>% 
#   summarise(n_outlier = n()) %>% 
#   view()

# ASHN_MYV_FST_outliers = inner_join(ASHN_Fst_clean, 
#                                  MYV_Fst_clean, 
#                                  by = c('CHR', 
#                                         'SNP', 
#                                         'POS'))
# 
# ASHN_SKR_FST_outliers = inner_join(ASHN_Fst_clean, 
#                                    SKR_Fst_clean, 
#                                    by = c('CHR', 
#                                           'SNP', 
#                                           'POS'))
# 
# MYV_SKR_FST_outliers = inner_join(MYV_Fst_clean, 
#                                    SKR_Fst_clean, 
#                                    by = c('CHR', 
#                                           'SNP', 
#                                           'POS'))
# 
# 
# ASHN_GTSCSWY_FST_outliers = inner_join(ASHN_Fst_clean, 
#                                    GTS_CSWY_Fst_clean, 
#                                    by = c('CHR', 
#                                           'SNP', 
#                                           'POS'))
# MYV_GTSCSWY_FST_outliers = inner_join(MYV_Fst_clean, 
#                                        GTS_CSWY_Fst_clean, 
#                                        by = c('CHR', 
#                                               'SNP', 
#                                               'POS'))
# SKR_GTSCSWY_FST_outliers = inner_join(SKR_Fst_clean, 
#                                        GTS_CSWY_Fst_clean, 
#                                        by = c('CHR', 
#                                               'SNP', 
#                                               'POS'))
# 

# common_FST_outliers = inner_join(ASHN_Fst_clean,
#            MYV_Fst_clean, 
#            by = c('CHR', 
#                   'SNP', 
#                   'POS')) %>% 
#   inner_join(., 
#              SKR_Fst_clean, 
#              by = c('CHR', 
#                     'SNP', 
#                     'POS')) %>% 
#   inner_join(., 
#              GTS_CSWY_Fst_clean, 
#              by = c('CHR', 
#                     'SNP', 
#                     'POS'))
# 
# 
# common_FST_outliers %>% 
#   select(CHR, 
#          SNP, 
#          POS, 
#          starts_with('FST_zero'))

# Fst sliding window ------------------------------------------------------

## per population FST outliers
# ASHN_Fst_clean = read_csv('ASHN_Fst_clean.csv') %>%
#   stickle_CHR_reorder() %>%
# dist_cal()
# MYV_Fst_clean = read_csv('MYV_Fst_clean.csv') %>%
#   stickle_CHR_reorder() %>%
#   dist_cal()
# SKR_Fst_clean = read_csv('SKR_Fst_clean.csv') %>%
#   stickle_CHR_reorder() %>%
#   dist_cal()
# GTS_CSWY_Fst_clean = read_csv('GTS_CSWY_Fst_clean.csv') %>%
#   stickle_CHR_reorder() %>%
#   dist_cal()

##Common FST outliers
WC_Fst_clean = read_csv('WC_Fst_clean.csv') %>%
  stickle_CHR_reorder() %>%
  dist_cal()

## 25kb
fst_25_position = winScan(x = WC_Fst_clean, 
                       groups = 'CHR', 
                       position = 'POS',
                       values = 'FST', 
                       win_size = 25000, 
                       win_step = 24999, 
                       funs = c('mean', 'sd'))

fst_25_position = fst_25_position %>%
  as_tibble() %>% 
  filter(FST_n >= 3) %>% 
  write_tsv('WC_Fst_25kb_3obs_window.txt')
## Write the txt file for each window size. 
## Need to compare the different window sizes to see which one
## is the most appropriate. 
## small window size == greater chance for false positives
## large window size == less chance to find differences

# write_tsv(fst_25_position, 
#           'MYV_Fst_25Kb_3obs_window.txt')
# 

## sliding window analysis for 50Kb windows
## with 1Kb overlap between windows
# fst_position = winScan(x = GTS_CSWY_Fst_clean, 
#                        groups = 'CHR', 
#                        position = 'POS',
#                        values = 'FST', 
#                        win_size = 50000, 
#                        win_step = 49000, 
#                        funs = c('mean', 'sd'))
# 
# fst_position = fst_position %>%
#   as_tibble() %>% 
#   filter(FST_n >= 3)
# ## Write the txt file for each window size. 
# ## Need to compare the different window sizes to see which one
# ## is the most appropriate. 
# ## small window size == greater chance for false positives
# ## large window size == less chance to find differences
# 
# write_tsv(fst_position, 
#           'GTS_CSWY_Fst_50Kb_3obs_window.txt')

##


# FST outliers 50Kb regions ------------------------------------------------
## already filtered for the FST_n containing at least 3 observations
# WC_50kb = read_tsv('WC_Fst_50Kb_3obs_window.txt')
# ASHN_50kb = read_tsv('ASHN_Fst_50Kb_3obs_window.txt')
# MYV_50kb = read_tsv('MYV_Fst_50Kb_3obs_window.txt')
# SKR_50kb = read_tsv('SKR_Fst_50Kb_3obs_window.txt')
# GTS_CSWY_50kb = read_tsv('GTS_CSWY_Fst_50Kb_3obs_window.txt')
WC_25kb = read_tsv('WC_Fst_25Kb_3obs_window.txt')
ASHN_25kb = read_tsv('ASHN_Fst_25Kb_3obs_window.txt')
MYV_25kb = read_tsv('MYV_Fst_25Kb_3obs_window.txt')
SKR_25kb = read_tsv('SKR_Fst_25Kb_3obs_window.txt')
GTS_CSWY_25kb = read_tsv('GTS_CSWY_Fst_25Kb_3obs_window.txt')

WC_25kb %>%
  SW_top_0.5_outliers() %>%
  # write_csv('SKR_25Kb_0.5%_Fst_outlier.csv')
  write_tsv('WC_25Kb_0.5%_Fst_outlier.tsv')

# ASHN_25kb %>%
#   SW_top_5_outliers() %>%
#   write_csv('GTS_CSWY_25Kb_Fst_outlier.csv')



# Fst 25kb outlier overlap ------------------------------------------------
WC_25_top5 = read_csv('WC_25Kb_0.5%_Fst_outlier.csv') 
MYV_25_top5 = read_csv('MYV_25Kb_0.5%_Fst_outlier.csv') 
ASHN_25_top5 = read_csv('ASHN_25Kb_0.5%_Fst_outlier.csv') 
SKR_25_top5 = read_csv('SKR_25Kb_0.5%_Fst_outlier.csv') 
GTS_CSWY_25_top5 = read_csv('GTS_CSWY_25Kb_0.5%_Fst_outlier.csv') 

intersect(ASHN_25_top5,
          MYV_25_top5, 
          by = c('CHR',
                 'win_mid'))

inner_join(GTS_CSWY_25_top5,
          WC_25_top5, 
          by = c('CHR',
                 'win_mid')) 

# FST 25kb manhattan plot -------------------------------------------------


# WC_25_top5 %>% 
#   filter(FST_mean >= 0.04) %>%
#   # filter(FST_mean >= 0.03) %>% 
#   # filter(FST_mean >= 0.018, 
#   #        CHR == 'chr_XIX') %>% 
#   View()
# location_cols = c('#00798c',
#                   '#003d5b',
#                   '#edae49',
#                   '#d1495b',
#                   '#30638e')

WC_25_top5 = read_csv('WC_25Kb_0.5%_Fst_outlier.csv') 
WC_25kb = read_tsv('WC_Fst_25Kb_3obs_window.txt') 

WC_25_window = Fst_manhatan_format(Fst_data = WC_25kb, 
                                   Fst_outliers = WC_25_top5) %>% 
  stickle_CHR_reorder() %>% 
  SW_dist_cal()
WC_25_axis_df = axis_df(WC_25_window)

outs = WC_25_window %>% 
  filter(value == 'Outlier')
neutral = WC_25_window %>% 
  filter(value == 'Neutral')

WC_25_Zoomed_plot = ggplot(neutral, 
       aes(x = BPcum, 
           y = FST_mean))+
  # plot the non outliers in grey
  geom_point(aes(color = as.factor(CHR)), 
             alpha = 0.8, 
             size = 1.3)+
  ## alternate colors per chromosome
  scale_color_manual(values = rep(c("grey", "dimgrey"), 39))+
  ## plot the outliers on top of everything
  ## currently digging this hot pink colour
  geom_point(data = outs,
             col = '#fb6f92',
             alpha=0.8, 
             size=1.3)+
  scale_x_continuous(label = WC_25_axis_df$CHR, 
                     breaks = WC_25_axis_df$center)+
  scale_y_continuous(expand = c(0, 0), 
                     limits = c(0,0.10))+
  labs(x = 'Cumulative base pair', 
       y = 'Fst', 
       title = 'B)')+
  theme(legend.position="none",
        # panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(), 
        axis.text.x = element_text(size = 9, 
                                   angle = 90), 
        axis.title = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12))

ggsave(file = 'Zoomed_common_adapt_stickleback_0.5%_FST_25KB_manhattan_plot.tiff', 
       path = '~/Parsons_Postdoc/Stickleback_Genomic/Figures/', 
       plot = WC_25_Zoomed_plot, 
       dpi = 'retina', 
       units = 'cm', 
       width = 30, 
       height = 15)

WC_25_region_man = Fst_manhattan(non_outs = neutral, 
                              outs = outs, 
                              axisdf = WC_25_axis_df, 
                              xval = BPcum, 
                              yval = FST_mean, 
                              chr = neutral$CHR,
                              out_col = '#fb6f92', 
                              plot_letter = 'A)')
 
ASHN_25_top5 = read_csv('ASHN_25Kb_0.5%_Fst_outlier.csv') 
ASHN_25kb = read_tsv('ASHN_Fst_25Kb_3obs_window.txt') 

ASHN_25_window = Fst_manhatan_format(Fst_data = ASHN_25kb, 
                                   Fst_outliers = ASHN_25_top5) %>% 
  stickle_CHR_reorder() %>% 
  SW_dist_cal()
ASHN_25_axis_df = axis_df(ASHN_25_window)

outs = ASHN_25_window %>% 
  filter(value == 'Outlier')
neutral = ASHN_25_window %>% 
  filter(value == 'Neutral')


ASHN_25_region_man = Fst_manhattan(non_outs = neutral, 
                                 outs = outs, 
                                 axisdf = ASHN_25_axis_df, 
                                 xval = BPcum, 
                                 yval = FST_mean, 
                                 chr = neutral$CHR,
                                 out_col = '#00798c', 
                                 plot_letter = 'A) Áshildarholtsvatn geothermal-ambient comparison')



MYV_25_top5 = read_csv('MYV_25Kb_0.5%_Fst_outlier.csv') 

MYV_25kb = read_tsv('MYV_Fst_25Kb_3obs_window.txt') 

MYV_25_window = Fst_manhatan_format(Fst_data = MYV_25kb, 
                                     Fst_outliers = MYV_25_top5) %>% 
  stickle_CHR_reorder() %>% 
  SW_dist_cal()
MYV_25_axis_df = axis_df(MYV_25_window)

outs = MYV_25_window %>% 
  filter(value == 'Outlier')
neutral = MYV_25_window %>% 
  filter(value == 'Neutral')


MYV_25_region_man = Fst_manhattan(non_outs = neutral, 
                                   outs = outs, 
                                   axisdf = MYV_25_axis_df, 
                                   xval = BPcum, 
                                   yval = FST_mean, 
                                   chr = neutral$CHR,
                                   out_col = '#d1495b', 
                                   plot_letter = 'B) Mývatn geothermal-ambient comparison')

SKR_25_top5 = read_csv('SKR_25Kb_0.5%_Fst_outlier.csv') 
SKR_25kb = read_tsv('SKR_Fst_25Kb_3obs_window.txt') 

SKR_25_window = Fst_manhatan_format(Fst_data = SKR_25kb, 
                                    Fst_outliers = SKR_25_top5) %>% 
  stickle_CHR_reorder() %>% 
  SW_dist_cal()
SKR_25_axis_df = axis_df(SKR_25_window)

outs = SKR_25_window %>% 
  filter(value == 'Outlier')
neutral = SKR_25_window %>% 
  filter(value == 'Neutral')


SKR_25_region_man = Fst_manhattan(non_outs = neutral, 
                                  outs = outs, 
                                  axisdf = SKR_25_axis_df, 
                                  xval = BPcum, 
                                  yval = FST_mean, 
                                  chr = neutral$CHR,
                                  out_col = '#30638e', 
                                  plot_letter = 'C) Sauðárkrókur geothermal-ambient comparison')


GTS_CSWY_25_top5 = read_csv('GTS_CSWY_25Kb_0.5%_Fst_outlier.csv') 
GTS_CSWY_25kb = read_tsv('GTS_CSWY_Fst_25Kb_3obs_window.txt') 

GTS_CSWY_25_window = Fst_manhatan_format(Fst_data = GTS_CSWY_25kb, 
                                    Fst_outliers = GTS_CSWY_25_top5) %>% 
  stickle_CHR_reorder() %>% 
  SW_dist_cal()
GTS_CSWY_25_axis_df = axis_df(GTS_CSWY_25_window)

outs = GTS_CSWY_25_window %>% 
  filter(value == 'Outlier')
neutral = GTS_CSWY_25_window %>% 
  filter(value == 'Neutral')


GTS_CSWY_25_region_man = Fst_manhattan(non_outs = neutral, 
                                  outs = outs, 
                                  axisdf = GTS_CSWY_25_axis_df, 
                                  xval = BPcum, 
                                  yval = FST_mean, 
                                  chr = neutral$CHR,
                                  out_col = '#edae49', 
                                  plot_letter = 'D) Grettislaug-Garðsvatn comparison')



Fst_region_combo = (ASHN_25_region_man|MYV_25_region_man)/(SKR_25_region_man|GTS_CSWY_25_region_man)|WC_25_region_man

## ggsave that plot

ggsave(file = 'stickleback_0.5%_FST_25KB_manhattan_plot.tiff', 
       path = '~/Parsons_Postdoc/Stickleback_Genomic/Figures/', 
       plot = Fst_region_combo, 
       dpi = 'retina', 
       units = 'cm', 
       width = 55, 
       height = 20)

local_adapt = (ASHN_25_region_man|MYV_25_region_man)/(SKR_25_region_man|GTS_CSWY_25_region_man)

ggsave(file = 'local_adapt_stickleback_0.5%_FST_25KB_manhattan_plot.tiff', 
       path = '~/Parsons_Postdoc/Stickleback_Genomic/Figures/', 
       plot = local_adapt, 
       dpi = 'retina', 
       units = 'cm', 
       width = 30, 
       height = 20)

common_adapt_zoom = WC_25_Zoomed_plot
ggsave(file = 'common_adapt_stickleback_0.5%_FST_25KB_manhattan_plot.tiff', 
       path = '~/Parsons_Postdoc/Stickleback_Genomic/Figures/', 
       plot = common_adapt_zoom, 
       dpi = 'retina', 
       units = 'cm', 
       width = 30, 
       height = 15)


common_adapt = WC_25_region_man
ggsave(file = 'common_adapt_stickleback_0.5%_FST_25KB_manhattan_plot.tiff', 
       path = '~/Parsons_Postdoc/Stickleback_Genomic/Figures/', 
       plot = common_adapt, 
       dpi = 'retina', 
       units = 'cm', 
       width = 30, 
       height = 15)


# WC_outliers_per_chr -----------------------------------------------------

WC_25_top5 = read_csv('WC_25Kb_0.5%_Fst_outlier.csv') 
WC_25kb = read_tsv('WC_Fst_25Kb_3obs_window.txt') 

WC_25_window = Fst_manhatan_format(Fst_data = WC_25kb, 
                                   Fst_outliers = WC_25_top5) %>% 
  stickle_CHR_reorder() %>% 
  SW_dist_cal()

outs = WC_25_window %>% 
  filter(value == 'Outlier')
neutral = WC_25_window %>% 
  filter(value == 'Neutral')


## Need the ratio of outlier loci to neutral loci across each chromsome

out_dat = outs %>%
  group_by(CHR) %>% 
  summarize(out_num_chr_out = n(), 
            mean_win_spot_out = mean(win_mid), 
            FST_mean_chr_out = mean(FST_mean))

neutral_dat = neutral %>% 
  group_by(CHR) %>% 
  summarize(out_num_chr_neutral = n(), 
            mean_win_spot_neutral = mean(win_mid), 
            FST_mean_chr_neutral = mean(FST_mean))

WC_25_window_sum = full_join(out_dat, 
          neutral_dat)

WC_25_window_sum %>% 
  group_by(CHR) %>% 
  mutate(ratio_outs_neutral = out_num_chr_out/out_num_chr_neutral) %>% 
  select(CHR, 
         FST_mean_chr_out, 
         FST_mean_chr_neutral, 
         ratio_outs_neutral) 
## Third most outliers compared to neutral snps

## But we want to see if it has the most SNPs within a given location
## That's what matters, are all of the outlier snps that we see within
## a specific genomic range along the chromosome. 

## We want the concentration of outliers along the chromosome. 

chr21_outs = outs %>% 
  filter(CHR == 'chr_XXI')

chr21_neut = neutral %>% 
  filter(CHR == 'chr_XXI')

chr21_axis_df = WC_25_axis_df %>% 
  filter(CHR == 'chr_XXI')

chr21_Fst_25kb_win = ggplot(chr21_neut, 
                           aes(x = win_mid, 
                               y = FST_mean))+
  geom_point(col = 'grey', 
             alpha = 0.8, 
             size = 2)+
  geom_point(data = chr21_outs,
             col = '#fb6f92',
             alpha=0.8, 
             size=2)+
  scale_y_continuous(expand = c(0, 0), 
                     limits = c(0,0.05))+
  labs(x = 'Base pair position', 
       y = 'Fst', 
       title = 'A)')+
  theme(legend.position="none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))

ggsave(file = 'Zoomed_chr21_0.5%_FST_25KB_manhattan_plot.tiff', 
       path = '~/Parsons_Postdoc/Stickleback_Genomic/Figures/', 
       plot = chr21_Fst_25kb_win, 
       dpi = 'retina', 
       units = 'cm', 
       width = 15, 
       height = 10)
##
# Chr XXI SNPS per population --------------------------------------------------

ASHN_Fst_clean = read_csv('ASHN_Fst_clean.csv') %>% 
  stickle_CHR_reorder() %>% 
  dist_cal()
MYV_Fst_clean = read_csv('MYV_Fst_clean.csv') %>% 
  stickle_CHR_reorder() %>% 
  dist_cal()
SKR_Fst_clean = read_csv('SKR_Fst_clean.csv') %>% 
  stickle_CHR_reorder() %>% 
  dist_cal()
GTS_CSWY_Fst_clean = read_csv('GTS_CSWY_Fst_clean.csv') %>% 
  stickle_CHR_reorder() %>% 
  dist_cal()
WC_Fst_clean = read_csv('WC_Fst_clean.csv') %>% 
  stickle_CHR_reorder() %>% 
  dist_cal()
## calculate the center of the chromosome
axisdf = axis_df(ASHN_Fst_clean)
axisdf = axis_df(MYV_Fst_clean)
axisdf = axis_df(SKR_Fst_clean)
axisdf = axis_df(GTS_CSWY_Fst_clean)
axisdf = axis_df(WC_Fst_clean)

non_outs = 
  WC_Fst_clean %>%
  # ASHN_Fst_clean %>%
  # MYV_Fst_clean %>%
  # SKR_Fst_clean %>%
  # GTS_CSWY_Fst_clean %>%
  filter(value == 'Neutral') %>% 
  filter(CHR == 'chr_XXI')

## Get the outliers
outs = 
  WC_Fst_clean %>%
  # ASHN_Fst_clean %>%
  # MYV_Fst_clean %>%
  # SKR_Fst_clean %>%
  # GTS_CSWY_Fst_clean %>%
  filter(value == 'Outlier') %>% 
  filter(CHR == 'chr_XXI')

## ASHN colour = #00798c
## MYV colour = #d1495b
## SKR colour = #30638e
## GTS_CSWY colour = #edae49
## WC colour = #439a86



ggplot(non_outs, 
       aes(x = POS, 
           y = FST_zero))+
  # plot the non outliers in grey
  geom_point(color = 'grey', 
             alpha = 0.8, 
             size = 1.3)+
  ## alternate colors per chromosome
  ## plot the outliers on top of everything
  ## currently digging this hot pink colour
  geom_point(data = outs,
             col = '#439a86',
             alpha=0.8, 
             size=1.3)+
  scale_y_continuous(expand = c(0, 0), 
                     limits = c(0,1.0))+
  # geom_hline(yintercept = 0.00043, 
  #            linetype = 2, 
  #            col = 'Black')+
  # ylim(0,1.0)+
  # scale_y_reverse(expand = c(0, 0))+
  # remove space between plot area and x axis
  labs(x = 'Base pair', 
       y = 'Fst', 
       title = 'E)')+
  theme(legend.position="none",
        # panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(), 
        axis.text.x = element_text(size = 9, 
                                   angle = 90), 
        axis.title = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12))

##
# FST 50kb region manhattan plot ------------------------------------------
WC_top5 = read_csv('WC_50Kb_Fst_outlier.csv') 
WC_50kb = read_tsv('WC_Fst_50Kb_3obs_window.txt') 

WC_window_df = Fst_manhatan_format(Fst_data = WC_50kb, 
                                   Fst_outliers = WC_top5) %>% 
  stickle_CHR_reorder() %>% 
  SW_dist_cal()
WC_axis_df = axis_df(WC_window_df)

outs = WC_window_df %>% 
  filter(value == 'Outlier')
neutral = WC_window_df %>% 
  filter(value == 'Neutral')


WC_region_man = Fst_manhattan(non_outs = neutral, 
                              outs = outs, 
                              axisdf = WC_axis_df, 
                              xval = BPcum, 
                              yval = FST_mean, 
                              chr = neutral$CHR,
                              out_col = '#ef233c', 
                              plot_letter = 'E)')


## FUCK there is drift in BP due to calculating the BPcum
## for two different datasets! Need to calculate this together and then
## split them based on the outlier and neutral labels. 
ASHN_50kb = read_tsv('ASHN_Fst_50Kb_3obs_window.txt') 
ASHN_top5 = read_csv('ASHN_50Kb_Fst_outlier.csv')
ASHN_window_df = Fst_manhatan_format(Fst_data = ASHN_50kb, 
                                     Fst_outliers = ASHN_top5) %>% 
  stickle_CHR_reorder() %>% 
  SW_dist_cal()
ASHN_axis_df = axis_df(ASHN_window_df)

outs = ASHN_window_df %>% 
  filter(value == 'Outlier')
neutral = ASHN_window_df %>% 
  filter(value == 'Neutral')

ASHN_region_man = Fst_manhattan(non_outs = neutral, 
                              outs = outs, 
                              axisdf = ASHN_axis_df, 
                              xval = BPcum, 
                              yval = FST_mean, 
                              chr = neutral$CHR,
                              out_col = '#06d6a0', 
                              plot_letter = 'A)')



MYV_50kb = read_tsv('MYV_Fst_50Kb_3obs_window.txt')
MYV_top5 = read_csv('MYV_50Kb_Fst_outlier.csv')
MYV_window_df = Fst_manhatan_format(Fst_data = MYV_50kb, 
                                    Fst_outliers = MYV_top5) %>% 
  stickle_CHR_reorder() %>% 
  SW_dist_cal()
MYV_axis_df = axis_df(MYV_window_df)
outs = MYV_window_df %>% 
  filter(value == 'Outlier')
neutral = MYV_window_df %>% 
  filter(value == 'Neutral')

MYV_region_man = Fst_manhattan(non_outs = neutral, 
                                outs = outs, 
                                axisdf = MYV_axis_df, 
                                xval = BPcum, 
                                yval = FST_mean, 
                                chr = neutral$CHR,
                                out_col = '#d62828', 
                                plot_letter = 'B)')


SKR_50kb = read_tsv('SKR_Fst_50Kb_3obs_window.txt') 
SKR_top5 = read_csv('SKR_50Kb_Fst_outlier.csv')
SKR_window_df = Fst_manhatan_format(Fst_data = SKR_50kb, 
                                    Fst_outliers = SKR_top5) %>% 
  stickle_CHR_reorder() %>% 
  SW_dist_cal()
SKR_axis_df = axis_df(SKR_window_df)
outs = SKR_window_df %>% 
  filter(value == 'Outlier')
neutral = SKR_window_df %>% 
  filter(value == 'Neutral')
SKR_region_man = Fst_manhattan(non_outs = neutral, 
                               outs = outs, 
                               axisdf = MYV_axis_df, 
                               xval = BPcum, 
                               yval = FST_mean, 
                               chr = neutral$CHR,
                               out_col = '#5f0f40', 
                               plot_letter = 'C)')


GTS_CSWY_top5 = read_csv('GTS_CSWY_50Kb_Fst_outlier.csv')
GTS_CSWY_50kb = read_tsv('GTS_CSWY_Fst_50Kb_3obs_window.txt') 
GTS_CSWY_window_df = Fst_manhatan_format(Fst_data = GTS_CSWY_50kb, 
                                         Fst_outliers = GTS_CSWY_top5) %>% 
  stickle_CHR_reorder() %>% 
  SW_dist_cal()

GTS_CSWY_axis_df = axis_df(GTS_CSWY_window_df)

outs = GTS_CSWY_window_df %>% 
  filter(value == 'Outlier')
neutral = GTS_CSWY_window_df %>% 
  filter(value == 'Neutral')

GTS_CSWY_region_man = Fst_manhattan(non_outs = neutral, 
                               outs = outs, 
                               axisdf = GTS_CSWY_axis_df, 
                               xval = BPcum, 
                               yval = FST_mean, 
                               chr = neutral$CHR,
                               out_col = '#264653', 
                               plot_letter = 'D)')


Fst_region_combo = (ASHN_region_man|MYV_region_man)/(SKR_region_man|GTS_CSWY_region_man)|WC_region_man


## ggsave that plot

ggsave(file = 'stickleback_FST_50KB_manhattan_plot.tiff', 
       path = 'C:/Stickleback_Genomic/Figures/', 
       plot = Fst_region_combo, 
       dpi = 'retina', 
       units = 'cm', 
       width = 40, 
       height = 20)


##
# LFMM Analysis -----------------------------------------------------------


temp_lfmm = lfmm( "stickleback_lfmm.lfmm", 
                "temp.env", 
                K = 5, 
                repetitions = 5, 
                project = "new")


# LFMM Manhattan plot -----------------------------------------------------

# map_test = read_tsv('stickleback_clean_ped.map', 
#                     col_names = F)
# ped_test = read_table2('stickleback_clean_ped.ped', 
#                     col_names = F)
# 
# convert = ped2geno('stickleback_clean_ped.ped', 
#                    'stickleback_data.geno')
# 
# 
# convert = geno2lfmm("stickleback_data.geno", 
#                     'stickleback_lfmm.lfmm')
# # 

LFMM_data = read_csv('~/Parsons_Postdoc/Stickleback_Genomic/lfmm/Stickleback_LFMM_temperature_qvalues.csv')

map_test = read_tsv('stickleback_clean_ped.map', 
                    col_names = c('CHR', 
                                  'SNP', 
                                  'Genetic_pos', 
                                  'POS')) 

LFMM_data = bind_cols(map_test, 
                      LFMM_data) %>% 
  rename(qvalue = value) %>% 
  stickle_CHR_reorder() %>% 
  dist_cal()

LFMM_axis_df = axis_df(LFMM_data)

non_outs = LFMM_data %>% 
  filter(qvalue >= 0.01) %>% 
  mutate(qvalue_trans = -log10(qvalue))

outs = LFMM_data %>% 
  filter(qvalue < 0.01) %>% 
  mutate(qvalue_trans = -log10(qvalue))

LFMM_Raw_qvalues = ggplot(non_outs, 
       aes(x = BPcum, 
           y = qvalue))+
  # plot the non outliers in grey
  geom_point(aes(color = as.factor(CHR)), 
             alpha = 0.8, 
             size = 1.3)+
  ## alternate colors per chromosome
  scale_color_manual(values = rep(c("grey", "dimgrey"), 39))+
  ## plot the outliers on top of everything
  ## currently digging this hot pink colour
  geom_point(data = outs,
             col = '#f72585',
             alpha=0.8, 
             size=1.3)+
  scale_x_continuous(label = LFMM_axis_df$CHR, 
                     breaks = LFMM_axis_df$center)+
  # scale_y_continuous(expand = c(0, 0))+
  # geom_hline(yintercept = 0.00043, 
  #            linetype = 2, 
  #            col = 'Black')+
  # ylim(0,1.0)+
  scale_y_reverse(expand = c(0, 0))+
  # remove space between plot area and x axis
  labs(x = 'Cumulative base pair', 
       y = 'Raw qvalues')+
  theme(legend.position="none",
        # panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(), 
        axis.text.x = element_text(size = 9, 
                                   angle = 90), 
        axis.title = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12))


LFMM_qvalue_logged = ggplot(non_outs, 
       aes(x = BPcum, 
           y = qvalue_trans))+
  # plot the non outliers in grey
  geom_point(aes(color = as.factor(CHR)), 
             alpha = 0.8, 
             size = 1.3)+
  ## alternate colors per chromosome
  scale_color_manual(values = rep(c("grey", "dimgrey"), 39))+
  ## plot the outliers on top of everything
  ## currently digging this hot pink colour
  geom_point(data = outs,
             col = '#f72585',
             alpha=0.8, 
             size=3)+
  scale_x_continuous(label = LFMM_axis_df$CHR, 
                     breaks = LFMM_axis_df$center)+
  scale_y_continuous(expand = c(0, 0))+
  # geom_hline(yintercept = 0.00043, 
  #            linetype = 2, 
  #            col = 'Black')+
  # ylim(0,1.0)+
  # scale_y_reverse(expand = c(0, 0))+
  # remove space between plot area and x axis
  labs(x = 'Cumulative base pair', 
       y = '-log(qvalue)')+
  theme(legend.position="none",
        # panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(), 
        axis.text.x = element_text(size = 9, 
                                   angle = 90), 
        axis.title = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12))

ggsave('~/Parsons_Postdoc/Stickleback_Genomic/Figures/LFMM_temp_log_qvalues.tiff', 
       plot = LFMM_qvalue_logged, 
       dpi = 'retina', 
       width = 10, 
       height = 5)

ggsave('~/Parsons_Postdoc/Stickleback_Genomic/Figures/LFMM_temp_raw_qvalues.tiff', 
       plot = LFMM_Raw_qvalues, 
       dpi = 'retina', 
       width = 20, 
       height = 10)


# LFMM FST COMBO ----------------------------------------------------------
LFMM_data = read_csv('~/Parsons_Postdoc/Stickleback_Genomic/lfmm/Stickleback_LFMM_temperature_qvalues.csv')

map_test = read_tsv('stickleback_clean_ped.map', 
                    col_names = c('CHR', 
                                  'SNP', 
                                  'Genetic_pos', 
                                  'POS')) 

LFMM_data = bind_cols(map_test, 
                      LFMM_data) %>% 
  rename(qvalue = value) %>% 
  stickle_CHR_reorder() %>% 
  dist_cal()

# write_csv(LFMM_data, 
#           'LFMM_Mapped_Full_Data.csv')

WC_Fst = read_tsv('Warm_Cold_Fst.fst') %>% 
  na.omit() %>% 
  mutate(FST_zero = if_else(FST < 0, 0, FST))

LFMM_FST_data = inner_join(WC_Fst, 
           LFMM_data, 
           by = c('CHR', 
                  'SNP', 
                  'POS'))

LFMM_FST_axis_df = axis_df(LFMM_FST_data)

non_outs = LFMM_FST_data %>% 
  filter(qvalue >= 0.01) %>% 
  mutate(qvalue_trans = -log10(qvalue))

outs = LFMM_FST_data %>% 
  filter(qvalue < 0.01) %>% 
  mutate(qvalue_trans = -log10(qvalue))


LFMM_FST = ggplot(non_outs, 
                            aes(x = BPcum, 
                                y = FST_zero))+
  # plot the non outliers in grey
  geom_point(aes(color = as.factor(CHR)), 
             alpha = 0.8, 
             size = 1.3)+
  ## alternate colors per chromosome
  scale_color_manual(values = rep(c("grey", "dimgrey"), 23))+
  ## plot the outliers on top of everything
  ## currently digging this hot pink colour
  geom_point(data = outs,
             col = '#f72585',
             alpha=0.8, 
             size=3)+
  scale_x_continuous(label = LFMM_FST_axis_df$CHR, 
                     breaks = LFMM_FST_axis_df$center)+
  scale_y_continuous(expand = c(0, 0), 
                     limits = c(0, 0.25))+
  # geom_hline(yintercept = 0.00043, 
  #            linetype = 2, 
  #            col = 'Black')+
  # ylim(0,1.0)+
  # scale_y_reverse(expand = c(0, 0))+
  # remove space between plot area and x axis
  labs(x = 'Cumulative base pair', 
       y = 'Fst')+
  theme(legend.position="none",
        # panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(), 
        axis.text.x = element_text(size = 9, 
                                   angle = 90), 
        axis.title = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12))


outs %>% 
  filter(FST_zero > 0.15) %>% 
  View()

ggsave('~/Parsons_Postdoc/Stickleback_Genomic/Figures/LFMM_temp_FST.tiff', 
       plot = LFMM_FST, 
       dpi = 'retina', 
       width = 10, 
       height = 5)

##
# af-vapeR ----------------------------------------------------------------

setwd('~/Parsons_Postdoc/Stickleback_Genomic/afvaper/')

# library(remotes)
# remotes::install_github("JimWhiting91/afvaper")
library(afvaper)


## need to create a popmap file
## which has individual id then the ecotype id
popmap = read_table2('stickle_filtered_chr_fix.ped', 
         col_names = F) %>% 
  dplyr::select(1:2)

# View(popmap)

popmap1 = popmap %>% 
  slice(1:16) %>% 
  separate(col = X1, 
           into = c('garbage', 
                    'Morph'), 
           sep = '-') %>% 
  separate(col = Morph, 
           into = c('Morph', 
                    'garbage'), 
           sep = '_') %>% 
  dplyr::select(Morph, 
                X2)

popmap2 = popmap %>% 
  slice(17:86) %>% 
  separate(col = X1, 
           into = c('Garbage', 
                    'Morph'), 
           sep = '-') %>% 
  separate(col = Morph, 
           into = c('Morph', 
                    'garbage'), 
           sep = '_') %>% 
  dplyr::select(Morph, 
                X2) 

popmap3 = popmap %>% 
  slice(87:109) %>% 
  separate(col = X1, 
           into = c('garbage', 
                    'Morph', 
                    'garbage2'), 
           sep = '_') %>% 
  dplyr::select(Morph, 
                X2)

popmap = bind_rows(popmap1, 
                   popmap2, 
                   popmap3) %>% 
  filter(Morph %in% c('MYVC', 
                      'MYVW', 
                      'ASHNC', 
                      'ASHNW', 
                      'SKRC', 
                      'SKRW'))

popmap = mutate(.data = popmap,
                      Ecotype = as.factor(case_when(
                        Morph == 'ASHNC' ~ 'Cold2',
                        Morph == 'ASHNW' ~ 'Warm2',
                        Morph == 'CSWY' ~ 'Cold4',
                        Morph == 'GTS' ~ 'Warm4',
                        Morph == 'MYVC' ~ 'Cold1',
                        Morph == 'MYVW' ~ 'Warm1',
                        Morph == 'SKRC' ~ 'Cold3',
                        Morph == 'SKRW' ~ 'Warm3'))) %>% 
  dplyr::select(X2, 
                Ecotype)

write_tsv(popmap, 
          'Stickleback_afvaper_round2_popmap.txt', 
          col_names = F)


fam_file = read_table2('strickle_No_GTS_CSWY_chr_fix.fam', 
                       col_names = F) %>% 
  dplyr::select(X1, X2)

popmap_file = read_tsv('Stickleback_afvaper_round2_popmap.txt', 
                       col_names = F)


popmap_file$X1==fam_file$X1

## Need to split the genomic data for each individual chromosome
## and then output that as a vcf file for use in plink

## Need to create a chromosome conversion function to change the stickleback
## custom chromosome names to straight numbers to be used in a forloop
## that way we always have an index of the what the chromosomes mean


# stickle_map = read_tsv('stickleback_maf0.05_ldpruned_filtered.map', 
#                        col_names = c('Chromosome', 
#                                      'SNP', 
#                                      'Genetic_pos', 
#                                      'Physical_pos'))

stickle_map = read_tsv('stickleback_No_GTS_CSWY.map', 
                       col_names = c('Chromosome', 
                                     'SNP', 
                                     'Genetic_pos', 
                                     'Physical_pos'))

Chr_convert(data = stickle_map) %>% 
  dplyr::select(chr_num, 
                SNP, 
                Genetic_pos, 
                Physical_pos) %>% 
  write_tsv('strickle_No_GTS_CSWY_chr_fix.map', 
            col_names = F)
  

## second we need to read in the fasta file to determine
## the number of null parameters to calculate per chromosome


chr_size = read_tsv('stickle_filtered_chr_fix.map', 
                    col_names = c('Chromosome', 
                                  'SNP', 
                                  'Genetic_pos', 
                                  'Physical_pos')) %>% 
# Chr_convert(data = stickle_map) %>% 
#   dplyr::select(chr_num, 
#                 SNP, 
#                 Genetic_pos, 
#                 Physical_pos) %>% 
  # stickle_map %>% 
  group_by(Chromosome) %>% 
  summarize(min_BP = min(Physical_pos), 
            max_BP = max(Physical_pos)) %>% 
  mutate(chr_size = max_BP - min_BP) %>% 
  # rename(Chromosome = chr_num) %>% 
  dplyr::select(Chromosome, 
                chr_size)

total_perms = 100000 
total_perms = 10000

chr_props = chr_size$chr_size/sum(chr_size$chr_size)
chr_perms = data.frame(chr = chr_size$Chromosome, 
                       perms = round(chr_props * total_perms))

# This gives us approximately 100000 null perms in total, 
##distributed across the genome according to relative size of chromosomes...
chr_perms %>% 
  as_tibble() %>% 
  arrange(chr) %>%
  write_tsv('Stickleback_Chromosome_sizes_100000.txt')






# afvaper results ---------------------------------------------------------

setwd('~/Parsons_Postdoc/Stickleback_Genomic/afvaper/afvaper results/')
# eig1_10snps = read_csv('afvaper_eigenvector1_results.csv')
# 
# View(eig1_10snps)
# 
# eig2_10snps = read_csv('afvaper_eigenvector2_results.csv')
# 
# View(eig2_10snps)


eig1_50snps = read_csv('afvaper_eigenvector1_results_50snp_window.csv')
# eig2_50snps = read_csv('afvaper_eigenvector2_results_50snp_window.csv')

View(eig1_50snps)
# View(eig2_50snps)

map = read_tsv('stickle_filtered_chr_fix.map', 
               col_names = c('CHR',
                             'SNP', 
                             'Genetic_pos', 
                             'Physical_pos'))

map %>% 
  # group_by(CHR) %>% 
  distinct(CHR, 
           .keep_all = T) %>% 
  View()

eig1_50snps = read_csv('afvaper_eigenvector1_results_50snp_window.csv')
# 
View(eig1_50snps)


setwd('~/Parsons_Postdoc/Stickleback_Genomic/vcf_filter/')

MYV_freq = read_tsv('MYV_afvaper_region_allele_freq.frq', 
                    col_names = c('CHR', 
                                  'SNP', 
                                  'A1', 
                                  'A2', 
                                  'Freq', 
                                  'Num_indivs'))




SKR_freq = read_tsv('SKR_afvaper_region_allele_freq.frq', 
                    col_names = c('CHR', 
                                  'SNP', 
                                  'A1', 
                                  'A2', 
                                  'Freq', 
                                  'Num_indivs'))

ASHN_freq = read_tsv('ASHN_afvaper_region_allele_freq.frq', 
                    col_names = c('CHR', 
                                  'SNP', 
                                  'A1', 
                                  'A2', 
                                  'Freq', 
                                  'Num_indivs'))

GTS_freq = read_tsv('GTS_afvaper_region_allele_freq.frq', 
                     col_names = c('CHR', 
                                   'SNP', 
                                   'A1', 
                                   'A2', 
                                   'Freq', 
                                   'Num_indivs'))

CSWY_freq = read_tsv('CSWY_afvaper_region_allele_freq.frq', 
                     col_names = c('CHR', 
                                   'SNP', 
                                   'A1', 
                                   'A2', 
                                   'Freq', 
                                   'Num_indivs'))

# 
# afvaper try chr1 --------------------------------------------------------

chr1_vcf = read.vcfR('stickle_filtered_1.vcf')

## individual ids got duplicated in the per chr vcfs somehow
chr1_vcf@gt

popmap = read_tsv('Stickleback_afvaper_popmap_dup.txt', 
                  col_names = F) %>% 
  as.data.frame()

# popmap = popmap %>%
#   mutate(X3 = X1) %>%
#   as_tibble() %>%
#   dplyr::select(X1,
#                 X3,
#                 X2)
# 
# popmap$X4 <- paste0(popmap$X1,
#                     "_",
#                     popmap$X3)
# 
# popmap %>%
#   select(X4,
#          X2) %>%
#   write_tsv('Stickleback_afvaper_popmap_dup.txt',
#             col_names = F)
# 

# unique(popmap[,2])

## Need to restructure the vector list into Warm1-cold1, 
## warm2-cold2 etc... with only a single pair it wont work, we 
## need to set each population as having a warm cold pairing for 
## the analysis to actually identify alleles that are parallel among
## all 4 axes/population pairs. 
# 
# input_vectors = list(pair1 = c("C1", "W1"),
#                      pair2 = c("C2", "W2"),
#                      pair3 = c("C3", "W3"),
#                      pair4 = c("C4", "W4"))

input_vectors = list(pair1 = c("Warm1", "Cold1"),
                     pair2 = c("Warm2", "Cold2"),
                     pair3 = c("Warm3", "Cold3"),
                     pair4 = c("Warm4", "Cold4"))


# Set our window size
window_snps = 50

# Calculate Allele Frequency Change Vector Matrices
AF_input = calc_AF_vectors(vcf = chr1_vcf,
                            window_size = window_snps,
                            popmap = popmap,
                            vectors = input_vectors,
                            n_cores = 1,
                            data_type = "vcf")




# How many permutations to run
null_perm_N = 1000

# Calculate Allele Frequency Change Vector Matrices
null_input = calc_AF_vectors(vcf = chr1_vcf,
                              window_size = window_snps,
                              popmap = popmap,
                              vectors = input_vectors,
                              n_cores = 1,
                              null_perms = null_perm_N,
                              data_type = "vcf")


chr_perms = read_tsv('Stickleback_Chromosome_sizes.txt')



# afvaper round2 results --------------------------------------------------

setwd('~/Parsons_Postdoc/Stickleback_Genomic/vcf_filter/')

afvaper_results = read_csv('afvaper_round2_eigenvector1_results_50snp_window.csv')
LFMM_data = read_csv('~/Parsons_Postdoc/Stickleback_Genomic/lfmm/Stickleback_LFMM_temperature_qvalues.csv')

map_test = read_tsv('stickleback_clean_ped.map', 
                    col_names = c('CHR', 
                                  'SNP', 
                                  'Genetic_pos', 
                                  'POS')) 

LFMM_data = bind_cols(map_test, 
                      LFMM_data) %>% 
  rename(qvalue = value) %>% 
  stickle_CHR_reorder() %>% 
  dist_cal()

# write_csv(LFMM_data, 
#           'LFMM_Mapped_Full_Data.csv')

WC_Fst = read_tsv('Warm_Cold_Fst.fst') %>% 
  na.omit() %>% 
  mutate(FST_zero = if_else(FST < 0, 0, FST))

LFMM_FST_data = inner_join(WC_Fst, 
                           LFMM_data, 
                           by = c('CHR', 
                                  'SNP', 
                                  'POS'))


LFMM_FST_data %>% 
  filter(CHR == 'chr_XXI', 
         POS >= 11044262, 
         POS <= 11213805) %>% 
  arrange(FST_zero) %>% 
  View()

LFMM_FST_data %>% 
  filter(CHR == 'chr_XXI', 
         POS >= 11044262, 
         POS <= 11213805) %>% 
  ggplot(aes(x = POS, 
             y = FST_zero))+
  geom_point()


LFMM_FST_data %>% 
  filter(CHR == 'chr_XXI', 
         POS >= 11370710, 
         POS <= 11574024) %>% 
  arrange(FST_zero) %>% 
  View()
LFMM_FST_data %>% 
  filter(CHR == 'chr_XXI', 
         POS >= 11370710, 
         POS <= 11574024) %>% 
  ggplot(aes(x = POS, 
             y = FST_zero))+
  geom_point()

# afvaper outlier overlap ----------------------------------------------------

afvaper_results = read_csv('afvaper_round2_eigenvector1_results_50snp_window.csv')

LFMM_data = read_csv('~/Parsons_Postdoc/Stickleback_Genomic/lfmm/Stickleback_LFMM_temperature_qvalues.csv')

map_test = read_tsv('stickleback_clean_ped.map', 
                    col_names = c('CHR', 
                                  'SNP', 
                                  'Genetic_pos', 
                                  'POS')) 

LFMM_data = bind_cols(map_test, 
                      LFMM_data) %>% 
  rename(qvalue = value) %>% 
  stickle_CHR_reorder() %>% 
  dist_cal()

# write_csv(LFMM_data, 
#           'LFMM_Mapped_Full_Data.csv')

WC_Fst = read_tsv('Warm_Cold_Fst.fst') %>% 
  na.omit() %>% 
  mutate(FST_zero = if_else(FST < 0, 0, FST))

pcadapt_outliers = read_csv('pcadapt_qvalues_all_snps.csv') %>% 
  rename(CHR = chromosome, 
         SNP = snp, 
         POS = physical_pos) %>% 
  inner_join(., 
            WC_Fst, 
            by = c('CHR', 
                   'SNP', 
                   'POS')) %>% 
  filter(qvalues <= 0.01) %>% 
  filter(CHR == 'chr_XXI')

LFMM_FST_data = inner_join(WC_Fst, 
                           LFMM_data, 
                           by = c('CHR', 
                                  'SNP', 
                                  'POS'))

WC_Fst_outs = read_csv('WC_Fst_clean.csv') %>% 
  stickle_CHR_reorder() %>% 
  dist_cal() %>% 
  filter(CHR == 'chr_XXI', 
         value == 'Outlier', 
         FST_zero >= 0.10)

LFMM_FST_data %>% 
  filter(CHR == 'chr_XXI', 
         POS >= 11044262, 
         POS <= 11574042) %>%  
  filter(qvalue <= 0.01)

LFMM_FST_data %>% 
  filter(CHR == 'chr_XXI', 
         POS >= 11044262, 
         POS <= 11574042) %>% 
  ggplot(aes(x = POS, 
             y = FST_zero))+
  geom_point()

LFMM_outs_xxi = LFMM_FST_data %>% 
filter(CHR == 'chr_XXI') %>% 
filter(qvalue <= 0.05)

LFMM_neutral_xxi = LFMM_FST_data %>% 
  filter(CHR == 'chr_XXI', 
         qvalue > 0.05)

ggplot(data = LFMM_neutral_xxi, 
       aes(x = POS, 
           y = FST_zero))+
  geom_point() +
  geom_point(data = WC_Fst_outs, 
             aes(x = POS, 
                 y = FST_zero), 
             col = '#3a86ff')+
  geom_point(data = LFMM_outs_xxi, 
             aes(x = POS, 
                 y = FST_zero), 
             col = '#f72585') +
  geom_point(data = pcadapt_outliers, 
             aes(x = POS, 
                 y = FST_zero), 
             col = '#6a994e')+
  geom_vline(xintercept = 11044262, 
             col = '#e63946', 
             linetype = 'dashed')+
  geom_vline(xintercept = 11213805, 
             col = '#e63946', 
             linetype = 'dashed')+
  geom_vline(xintercept = 11370710, 
             col = '#e63946', 
             linetype = 'dashed')+
  geom_vline(xintercept = 11574024, 
             col = '#e63946', 
             linetype = 'dashed') +
  labs(x = 'Physical position', 
       y = 'FST')+
  theme_bw()+
  theme(panel.grid = element_blank())

#
# Gene overlap ------------------------------------------------------------

setwd('~/Parsons_Postdoc/Stickleback_Genomic/')

afvaper_genes = read_tsv('afvapor_parallel_allele_freq_genes.txt') %>% 
  rename(Stickleback_Combo_Genes = 1)

all_other_genes = read_tsv('Stickleback_Combo_Genes.txt')

## Well shit, that didn't work at all
inner_join(afvaper_genes, 
           all_other_genes)
