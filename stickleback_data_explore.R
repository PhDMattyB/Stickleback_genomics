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
library(LEA)
library(reshape2)
library(qvalue)
library(tidyverse)
library(umap)
library(vcfR)
library(Rcpp)
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
                        population == 'ASHNC' ~ 'ASH',
                        population == 'ASHNW' ~ 'ASH',
                        population == 'CSWY' ~ 'CSWY',
                        population == 'GTS' ~ 'GTS',
                        population == 'MYVC' ~ 'MYV',
                        population == 'MYVW' ~ 'MYV',
                        population == 'SKRC' ~ 'SKR',
                        population == 'SKRW' ~ 'SKR')))



##
# pca plot ----------------------------------------------------------------

theme_set(theme_bw())

cold_warm_cols = c('#264653',
         '#e63946', 
         '#386641', 
         '#6a994e', 
         '#457b9d', 
         '#e76f51',
         '#6d6875', 
         '#c1121f')


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
  geom_point(aes(col = Location),
             size = 2)+
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

ggsave(file = 'stickleback_pca_cold_warm_split.tiff', 
       path = 'C:/Stickleback_Genomic/Figures/', 
       plot = stickleback_pca, 
       dpi = 'retina', 
       units = 'cm', 
       width = 20.0, 
       height = 13)


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
setwd('C:/Stickleback_Genomic/vcf_filter/')

convert = ped2geno('stickleback_maf0.05_ldpruned_filtered.ped',
                   'stickleback_data.geno')

# snmf('stickleback_data.geno',
#      K = 1:10,
#      entropy = T,
#      repetitions = 5,
#      project = 'new')

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
# head(qmatrix)
# dim(qmatrix)
## This will show the clusters each individual was assigned to
## The order is the same as the .ped file!!
## Take the first two columns of the ped file to group by population
apply(qmatrix, 1, which.max) %>%
  as_tibble() %>%
  dplyr::rename(Genetic_group = value) %>%
  write_tsv('stickleback_snmf_k5.txt')

## shitty base R plot
plot = barplot(qmatrix, 
               border = NA, 
               space = 0, 
               # col = my.colors, 
               # xlab = "Individuals",
               ylab = "Ancestry proportions", 
               main = "Ancestry matrix")

# bp = LEA::barchart(qmatrix, 
#                 K = 23, 
#                 border = NA, 
#                 space = 0, 
#                 xlab = 'Individuals', 
#                 ylab = 'Ancestry proportions', 
#                 main = 'Ancestry matrix')

as_tibble(qmatrix) %>%
  dplyr::rename(Q1 = V1,
                Q2 = V2,
                Q3 = V3,
                Q4 = V4, 
                Q5 = V5) %>% 
  write_csv('stickleback_snmf_qvalues_k5.csv')


snmf_data = read_csv('stickleback_snmf_qvalues_k5.csv')

identifiers

snmf_data = bind_cols(identifiers, 
                      snmf_data) %>% 
  # arrange(population, 
  #         individual_id)
  arrange(population)

snmf_melted = melt(snmf_data, 
                   id.vars = c('population', 
                               'individual_id')) %>% 
  as_tibble()

## need colour scheme
## real
location_cols = c('#264653',
                  '#5f0f40',
                  '#06d6a0',
                  '#219ebc',
                  '#d62828')

## snmf plot
snmf_plot = ggplot(data = snmf_melted, 
                   aes(x = individual_id,
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

WC_top_dist = WC_Fst[WC_Fst$FST_zero > quantile(WC_Fst$FST_zero, 
                                                      prob = 1-5/100),]

# ASHN_top_dist %>% write_csv('ASHN_FST_Outliers.csv')
WC_top_dist %>% 
  # group_by(CHR) %>% 
  summarise(min_fst = min(FST_zero), 
            max_fst = max(FST_zero), 
            mean_fst = mean(FST_zero)) 


## snps that are the top 5% fst distribution
ASHN_top_dist = ASHN_Fst[ASHN_Fst$FST_zero > quantile(ASHN_Fst$FST_zero, 
                                    prob = 1-5/100),]

# ASHN_top_dist %>% write_csv('ASHN_FST_Outliers.csv')
ASHN_top_dist %>% 
  # group_by(CHR) %>% 
  summarise(min_fst = min(FST_zero), 
            max_fst = max(FST_zero), 
            mean_fst = mean(FST_zero)) 

MYV_top_dist = MYV_Fst[MYV_Fst$FST_zero > quantile(MYV_Fst$FST_zero, 
                                                      prob = 1-5/100),]
# MYV_top_dist %>% write_csv('MYV_FST_Outliers.csv')
MYV_top_dist %>% 
  # group_by(CHR) %>% 
  summarise(min_fst = min(FST_zero), 
            max_fst = max(FST_zero), 
            mean_fst = mean(FST_zero))

SKR_top_dist = SKR_Fst[SKR_Fst$FST_zero > quantile(SKR_Fst$FST_zero, 
                                                   prob = 1-5/100),]
# SKR_top_dist %>% write_csv('SKR_FST_Outliers.csv')

SKR_top_dist %>% 
  # group_by(CHR) %>% 
  summarise(min_fst = min(FST_zero), 
            max_fst = max(FST_zero), 
            mean_fst = mean(FST_zero))

GTS_CSWY_top_dist = GTS_CSWY_Fst[GTS_CSWY_Fst$FST_zero > quantile(GTS_CSWY_Fst$FST_zero, 
                                                   prob = 1-5/100),]
# GTS_CSWY_top_dist %>% write_csv('GTS_CSWY_FST_Outliers.csv')

GTS_CSWY_top_dist %>% 
  # group_by(CHR) %>% 
  summarise(min_fst = min(FST_zero), 
            max_fst = max(FST_zero), 
            mean_fst = mean(FST_zero))


##
# FST distribution plots --------------------------------------------------

location_cols = c('#06d6a0',
                  '#264653',
                  '#219ebc',
                  '#d62828',
                  '#5f0f40')

ASHN_Fst_dist_plot = ASHN_Fst %>% 
  ggplot()+
  geom_density(aes(x = FST_zero), 
               col = '#06d6a0', 
               fill = '#06d6a0')+
  geom_density(data = ASHN_top_dist, 
               aes(x = FST_zero),
               col = '#000000',
               fill = '#000000')+
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
               col = '#d62828', 
               fill = '#d62828')+
  geom_density(data = MYV_top_dist, 
               aes(x = FST_zero),
               col = '#000000',
               fill = '#000000')+
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
               col = '#5f0f40', 
               fill = '#5f0f40')+
  geom_density(data = SKR_top_dist, 
               aes(x = FST_zero),
               col = '#000000',
               fill = '#000000')+
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
               col = '#264653', 
               fill = '#264653')+
  geom_density(data = GTS_CSWY_top_dist, 
               aes(x = FST_zero),
               col = '#000000',
               fill = '#000000')+
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
               col = '#ef233c', 
               fill = '#ef233c')+
  geom_density(data = WC_top_dist, 
               aes(x = FST_zero),
               col = '#000000',
               fill = '#000000')+
  labs(x = 'Fst', 
       y = 'Density', 
       title = 'E)')+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12))


Fst_dist_combo = (ASHN_Fst_dist_plot|MYV_Fst_dist_plot)/(SKR_Fst_dist_plot|GTS_CSWY_Fst_dist_plot) | WC_Fst_dist_plot 

ggsave(file = 'stickleback_FST_Distribution_plots.tiff', 
       path = 'C:/Stickleback_Genomic/Figures/', 
       plot = Fst_dist_combo, 
       dpi = 'retina', 
       units = 'cm', 
       width = 20.0, 
       height = 15)

# Fst distribution per chromosome plot ------------------------------------
location_cols = c('#06d6a0',
                  '#264653',
                  '#219ebc',
                  '#d62828',
                  '#5f0f40')

ASHN_per_chrom = ASHN_Fst %>% 
  ggplot()+
  geom_density(aes(x = FST_zero), 
               col = '#06d6a0', 
               fill = '#06d6a0')+
  geom_density(data = ASHN_top_dist, 
               aes(x = FST_zero),
               col = '#000000',
               fill = '#000000')+
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
               col = '#d62828', 
               fill = '#d62828')+
  geom_density(data = MYV_top_dist, 
               aes(x = FST_zero),
               col = '#000000',
               fill = '#000000')+
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
               col = '#5f0f40', 
               fill = '#5f0f40')+
  geom_density(data = SKR_top_dist, 
               aes(x = FST_zero),
               col = '#000000',
               fill = '#000000')+
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
               col = '#264653', 
               fill = '#264653')+
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
               col = '#ef233c', 
               fill = '#ef233c')+
  geom_density(data = WC_top_dist, 
               aes(x = FST_zero),
               col = '#000000',
               fill = '#000000')+
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



ggsave(file = 'stickleback_FST_Distribution_per_chrome.tiff', 
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
axisdf = axis_df(ASHN_Fst_clean)
axisdf = axis_df(MYV_Fst_clean)
axisdf = axis_df(SKR_Fst_clean)
axisdf = axis_df(GTS_CSWY_Fst_clean)
axisdf = axis_df(WC_Fst_clean)

non_outs = 
  # WC_Fst_clean %>% 
  # ASHN_Fst_clean %>%
  # MYV_Fst_clean %>%
  # SKR_Fst_clean %>%
  GTS_CSWY_Fst_clean %>%
  filter(value == 'Neutral') 

## Get the outliers
outs = 
  # WC_Fst_clean %>% 
  # ASHN_Fst_clean %>%
  # MYV_Fst_clean %>%
  # SKR_Fst_clean %>%
  GTS_CSWY_Fst_clean %>%
  filter(value == 'Outlier') 


## ASHN colour = #06d6a0
## MYV colour = #d62828
## SKR colour = #5f0f40
## GTS_CSWY colour = #264653
## WC colour = #ef233c

WC_Fst_manhattan = Fst_manhattan(non_outs = non_outs, 
                                 outs = outs, 
                                 axisdf = axisdf, 
                                 xval = BPcum, 
                                 yval = FST_zero, 
                                 chr = non_outs$CHR, 
                                 out_col = '#ef233c', 
                                 plot_letter = 'E)')


ASHN_Fst_manhattan = Fst_manhattan(non_outs = non_outs, 
                                       outs = outs, 
                                       axisdf = axisdf, 
                                       xval = BPcum, 
                                       yval = FST_zero, 
                                       chr = non_outs$CHR, 
                                       out_col = '#06d6a0', 
                                       plot_letter = 'A)')

MYV_Fst_manhattan = Fst_manhattan(non_outs = non_outs, 
                                   outs = outs, 
                                   axisdf = axisdf, 
                                   xval = BPcum, 
                                   yval = FST_zero, 
                                   chr = non_outs$CHR, 
                                   out_col = '#d62828', 
                                   plot_letter = 'B)')

SKR_Fst_manhattan = Fst_manhattan(non_outs = non_outs, 
                                  outs = outs, 
                                  axisdf = axisdf, 
                                  xval = BPcum, 
                                  yval = FST_zero, 
                                  chr = non_outs$CHR, 
                                  out_col = '#5f0f40', 
                                  plot_letter = 'C)')

GTS_CSWY_Fst_manhattan = Fst_manhattan(non_outs = non_outs,
              outs = outs, 
              axisdf = axisdf, 
              xval = BPcum, 
              yval = FST_zero, 
              chr = non_outs$CHR, 
              out_col = '#264653', 
              plot_letter = 'D)')


Fst_man_combo = (ASHN_Fst_manhattan|MYV_Fst_manhattan)/(SKR_Fst_manhattan|GTS_CSWY_Fst_manhattan)|WC_Fst_manhattan


## ggsave that plot

ggsave(file = 'stickleback_FST_manhattan_plot.tiff', 
       path = 'C:/Stickleback_Genomic/Figures/', 
       plot = Fst_man_combo, 
       dpi = 'retina', 
       units = 'cm', 
       width = 40, 
       height = 20)




# Outlier overlap ---------------------------------------------------------

## per population FST outliers
# ASHN_Fst_clean = read_csv('ASHN_Fst_clean.csv') %>% 
#   stickle_CHR_reorder() %>% 
#   dist_cal() %>% 
#   filter(value == 'Outlier')
# MYV_Fst_clean = read_csv('MYV_Fst_clean.csv') %>% 
#   stickle_CHR_reorder() %>% 
#   dist_cal()%>% 
#   filter(value == 'Outlier')
# SKR_Fst_clean = read_csv('SKR_Fst_clean.csv') %>% 
#   stickle_CHR_reorder() %>% 
#   dist_cal()%>% 
#   filter(value == 'Outlier')
# GTS_CSWY_Fst_clean = read_csv('GTS_CSWY_Fst_clean.csv') %>% 
#   stickle_CHR_reorder() %>% 
#   dist_cal()%>% 
#   filter(value == 'Outlier')

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
  filter(qvalue <= 0.01) %>% 
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
GTS_CSWY_Fst_clean = read_csv('GTS_CSWY_Fst_clean.csv') %>%
  stickle_CHR_reorder() %>%
  dist_cal()

##Common FST outliers
# WC_Fst_clean_all = read_csv('WC_Fst_clean.csv') %>%
#   stickle_CHR_reorder() %>%
#   dist_cal()

## 25kb
fst_25_position = winScan(x = GTS_CSWY_Fst_clean, 
                       groups = 'CHR', 
                       position = 'POS',
                       values = 'FST', 
                       win_size = 25000, 
                       win_step = 24999, 
                       funs = c('mean', 'sd'))

fst_25_position = fst_25_position %>%
  as_tibble() %>% 
  filter(FST_n >= 3) %>% 
  write_tsv('GTS_CSWY_Fst_25kb_3obs_window.txt')
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
# WC_25kb = read_tsv('WC_Fst_25Kb_3obs_window.txt')
ASHN_25kb = read_tsv('ASHN_Fst_25Kb_3obs_window.txt')
MYV_25kb = read_tsv('MYV_Fst_25Kb_3obs_window.txt')
SKR_25kb = read_tsv('SKR_Fst_25Kb_3obs_window.txt')
GTS_CSWY_25kb = read_tsv('GTS_CSWY_Fst_25Kb_3obs_window.txt')
GTS_CSWY_25kb %>% 
  SW_top_5_outliers() %>% 
  write_csv('GTS_CSWY_25Kb_Fst_outlier.csv')



# Fst 25kb outlier overlap ------------------------------------------------
WC_25_top5 = read_csv('WC_25Kb_Fst_outlier.csv') 
MYV_25_top5 = read_csv('MYV_25Kb_Fst_outlier.csv') 
ASHN_25_top5 = read_csv('ASHN_25Kb_Fst_outlier.csv') 
SKR_25_top5 = read_csv('SKR_25Kb_Fst_outlier.csv') 
GTS_CSWY_25_top5 = read_csv('GTS_CSWY_25Kb_Fst_outlier.csv') 

intersect(GTS_CSWY_25_top5,
          SKR_25_top5, 
          by = c('CHR',
                 'win_mid'))


# FST 25kb manhattan plot -------------------------------------------------


WC_25_top5 %>% 
  filter(FST_mean >= 0.04) %>%
  # filter(FST_mean >= 0.03) %>% 
  # filter(FST_mean >= 0.018, 
  #        CHR == 'chr_XIX') %>% 
  View()

WC_25_top5 = read_csv('WC_25Kb_Fst_outlier.csv') 
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


WC_25_region_man = Fst_manhattan(non_outs = neutral, 
                              outs = outs, 
                              axisdf = WC_25_axis_df, 
                              xval = BPcum, 
                              yval = FST_mean, 
                              chr = neutral$CHR,
                              out_col = '#ef233c', 
                              plot_letter = 'E)')

ASHN_25_top5 = read_csv('ASHN_25Kb_Fst_outlier.csv') 
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
                                 out_col = '#06d6a0', 
                                 plot_letter = 'A)')



MYV_25_top5 = read_csv('MYV_25Kb_Fst_outlier.csv') 

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
                                   out_col = '#d62828', 
                                   plot_letter = 'B)')

SKR_25_top5 = read_csv('SKR_25Kb_Fst_outlier.csv') 
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
                                  out_col = '#5f0f40', 
                                  plot_letter = 'C)')


GTS_CSWY_25_top5 = read_csv('GTS_CSWY_25Kb_Fst_outlier.csv') 
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
                                  out_col = '#264653', 
                                  plot_letter = 'D)')



Fst_region_combo = (ASHN_25_region_man|MYV_25_region_man)/(SKR_25_region_man|GTS_CSWY_25_region_man)|WC_25_region_man


## ggsave that plot

ggsave(file = 'stickleback_FST_25KB_manhattan_plot.tiff', 
       path = 'C:/Stickleback_Genomic/Figures/', 
       plot = Fst_region_combo, 
       dpi = 'retina', 
       units = 'cm', 
       width = 40, 
       height = 20)


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
                        # Morph == 'CSWY' ~ 'Cold4',
                        # Morph == 'GTS' ~ 'Warm4',
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

# eig1_50snps = read_csv('afvaper_eigenvector1_results_50snp_window.csv')
# 
# View(eig1_50snps)
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

# afvaper lfmm overlap ----------------------------------------------------

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
         POS <= 11574042) %>% View() 
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
