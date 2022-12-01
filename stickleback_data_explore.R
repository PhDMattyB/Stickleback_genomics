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
library(qvalue)
library(tidyverse)
library(umap)

# 


# pcadapt analysis --------------------------------------------------------

setwd('C:/Stickleback_Genomic/vcf_filter/')

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
  # geom_point(aes(col = Location),
  #            size = 2)+
  geom_point(aes(col = population),
             size = 2)+
  scale_color_manual(values = cold_warm_cols)+
  # scale_color_manual(values = location_cols)+
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



# snmf analysis -----------------------------------------------------------
setwd('C:/Stickleback_Genomic/vcf_filter/')

library(LEA)
library(viridis)
library(reshape2)

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
# rda analysis ------------------------------------------------------------


# Fst Plink ---------------------------------------------------------------

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

## Need to split based on each comparison 
## ASHW vs ASHC
## MYVW vs MYVC
## SKRW vs SKRC
## GTS vs CSWY

identifiers %>% 
  filter(population %in% c('ASHNW', 
                         'ASHNC')) %>%
  rename(`#population` = population) %>% 
  write_tsv('ASHN_Fst_grouping.txt')

## Holy fuck!! Make sure to use the actual family and individual
## identifiers in the fucking ped file. WOW

ped_ids = read_table2('stickleback_maf0.05_ldpruned_filtered.fam', 
                     col_names = F) %>%
              dplyr::select(X1,
                            X2) 

## Need to make a ped and map file for each of these comparisons
## the ped file is waaaay to big to open in R
## use the --keep or --keep-fam flags in plink to filter the 
## populations out. 
## the --keep file needs to be a text file with family and individual
## identifiers

ped_ids = bind_cols(ped_ids, 
          identifiers)


ped_ids %>% 
  filter(population %in% c('GTS', 
                           'CSWY')) %>%
  # rename(`#population` = population) %>%
  select(1:2) %>% 
  write_tsv('GTS_CSWY_Fst_keep.txt', 
            col_names = F)
