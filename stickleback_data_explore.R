##############################
## stickleback genomics explore
##
## Matt Brachmann (PhDMattyB)
##
## 2022-11-23
##
##############################

# library(patchwork)
# library(janitor)
# library(devtools)
# library(skimr)
# library(rsed)
# library(data.table)
# library(sjPlot)
library(tidyverse)
# library(vcfR)
library(pcadapt)
library(qvalue)
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

identifiers = read_table2('stickleback_maf0.05-ldpruned_nomissing.fam', 
         col_names = F) %>% 
  dplyr::select(X1) %>% 
  rename(population = X1) %>% 
  separate(col = population, 
           into = c('individual_id', 
                'population'), 
           sep = '-') %>% 
  separate(col = individual_id, 
            into = c('garbage', 
                   'temp_pop'), 
            sep = 'Sample_')
## going to need to separate into two data frames and 
## then bring them back together due to sample name issues

id_one = identifiers %>% 
  slice(1:86) %>% 
  select(population) %>% 
  separate(col = population, 
           into = c('population', 
                    'garbage'), 
           sep = '_') %>% 
  select(population)

id_two = identifiers %>% 
  slice(87:109) %>% 
  select(temp_pop) %>% 
  separate(col = temp_pop, 
           into = c('population', 
                    'garbage1', 
                    'garbage2'), 
           sep = '_') %>% 
  select(population)

identifiers = bind_rows(id_one, 
                        id_two)

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

ggsave(file = 'stickleback_pca.tiff', 
       path = 'C:/Stickleback_Genomic/Figures/', 
       plot = stickleback_pca, 
       dpi = 'retina', 
       units = 'cm', 
       width = 20.0, 
       height = 13)



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

# final_df  %>%  mutate(logpval = -log(pvalues), 
#                       logqval = -log(qvalues))

## read in the data
# final_df = read_csv('~/Charr_Adaptive_Introgression/Charr_Project_1/GeneticData/Charr_Lab_PCARDA_FinalDF.csv')

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

## get the absolute score for the RDA outliers
## the negatives plot like shit
# dist_cal$RDA_score_abs = abs(dist_cal$RDA_score)

# write_csv(dist_cal,
#           '~/Charr_Adaptive_Introgression/Charr_Project_1/RDA_Environmental_Variables/RDA_Outliers_distcal_df.csv')
# write_csv(axisdf,
#           '~/Charr_Adaptive_Introgression/Charr_Project_1/RDA_Environmental_Variables/RDA_Outliers_axisdf_df.csv')

## Get the neutral snps
non_outs = dist_cal %>% 
  filter(qvalues >= 0.05) %>% 
  mutate
## Get the outliers
outs = dist_cal %>% 
  filter(qvalues < 0.05)

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

plot(project,
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
  write_tsv('stickleback_k5.txt')

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
  write_csv('stickleback_snmf_qvalues_k4.csv')


snmf_data = read_csv('stickleback_snmf_qvalues_k4.csv')

snmf_data %>% 
  group_by(`#Familyid`) %>% 
  summarise(n = n())


theme_set(theme_bw())

library(viridis)
## plot is currently arranged by Latitude 
## See arrange function in snmf_data

## snmf latitude plot
snmf_plot = ggplot(data = snmf_data, 
                   aes(x = reorder(Individualid, Lat),
                       y = value, 
                       fill = variable, 
                       group = Lat))+
  geom_bar(stat = "identity", 
           width = 1)+
  scale_fill_manual(values = magma(n = 20))+
  labs(x = 'Individuals', 
       y = 'Ancestry proportion')+
  theme(axis.text.y = element_text(color = 'black'),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        # axis.text.x = element_text(angle = 90,
        #                            hjust = 1,
        #                            vjust = -0.09,
        #                            size = 6,
        #                            color = 'black'),
        legend.position = 'none')+
  scale_y_continuous(expand = c(0,0))+
  scale_x_discrete(guide = guide_axis(n.dodge = 5))

snmf_plot

ggsave('Mito_Nuclear_snmf_k20.tiff',
       # path = '~/BradburyLab_Postdoc/Charr_Project_1/Figures',
       plot = snmf_plot, 
       dpi = 'retina', 
       units = 'cm')

# rda analysis ------------------------------------------------------------

