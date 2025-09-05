##############################
## allelic heatmap chr 21
##
## Matt Brachmann (PhDMattyB)
##
## 26.08.2025
##
##############################


setwd('~/Parsons_Postdoc/Stickleback_Genomic/vcf_filter/')

library(tidyverse)
library(vcfR)
# library(ChromHeatMap)

vcf = read.vcfR('stickleback_filtered_vcf.vcf')

head(vcf)

chrxxi_data = vcf[getCHROM(vcf) == 'chr_XXI']

write.vcf(x = chrxxi_data, 
          'chrxxi_vcf.vcf')

head(chrxxi_data)

chrxxi_data@gt
chrxxi_data@fix

chr_xxi_meta_data = chrxxi_data@fix %>% 
  as.data.frame() %>% 
  as_tibble() %>% 
  dplyr::select(CHROM,
                POS, 
                ID, 
                REF, 
                ALT)

chr_xxi_geno_data = chrxxi_data@gt %>% 
  as.data.frame() %>% 
  as_tibble() %>% 
  dplyr::select(-FORMAT)


chr_XXI_SNPs = bind_cols(chr_xxi_meta_data, 
          chr_xxi_geno_data)



View(chr_XXI_SNPs)




# braf gene snp changes ---------------------------------------------------
vcf = read.vcfR('stickleback_filtered_vcf.vcf')


genome_meta_data = vcf@fix %>% 
  as.data.frame() %>% 
  as_tibble() %>% 
  dplyr::select(CHROM,
                POS, 
                ID, 
                REF, 
                ALT)


genome_meta_data %>% 
  filter(ID == 'chr_IV_22813388')
## A:T in genome browser

genome_meta_data %>% 
  filter(ID == 'chr_IV_22817061')
## T:A in genome browser

genome_meta_data %>% 
  filter(ID == 'chr_IV_22821293')
## T:A in genome browser

# pca on the chr xxi ------------------------------------------------------
library(vegan)
library(dartRverse)

chrxxi_dart = dartR.base::gl.read.vcf('chrxxi_vcf.vcf')

chrxxi_dart@other$loc.metrics$position = chrxxi_dart@position


filtered_gl_data = dartR.base::gl.filter.locmetric(
  x = chrxxi_dart,
  metric = "position", # The name of the metric you added
  lower = 9963830,
  upper = 11574024,
  keep = "within",
  verbose = 3 # For detailed output
)

pca_inversion = adegenet::glPca(filtered_gl_data, 
                      nf = 10)

pca_inversion$loadings
pca_inversion$scores
pca_inversion$eig

inversion_pca_load = pca_inversion$loadings %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  as.tibble() %>% 
  rename(SNP = rowname)

inversion_pca_score = pca_inversion$scores %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  as.tibble() %>% 
  rename(Individual = rowname) 


metadata = read_csv('stickleback_identifiers.csv')
pca_data = bind_cols(metadata, 
                     inversion_pca_score) %>% 
  mutate(.data = .,
                  Ecotype = as.factor(case_when(
                    population == 'ASHNC' ~ 'Cold',
                    population == 'ASHNW' ~ 'Warm',
                    population == 'CSWY' ~ 'Cold',
                    population == 'GTS' ~ 'Warm',
                    population == 'MYVC' ~ 'Cold',
                    population == 'MYVW' ~ 'Warm',
                    population == 'SKRC' ~ 'Cold',
                    population == 'SKRW' ~ 'Warm')))%>% 
  mutate(.data = .,
         POP_PAIR = as.factor(case_when(
           population == 'ASHNC' ~ 'ASHN',
           population == 'ASHNW' ~ 'ASHN',
           population == 'CSWY' ~ 'GTS-GAR',
           population == 'GTS' ~ 'GTS-GAR',
           population == 'MYVC' ~ 'MYV',
           population == 'MYVW' ~ 'MYV',
           population == 'SKRC' ~ 'SKR',
           population == 'SKRW' ~ 'SKR')))

plot_pal = c('#003049', 
             '#c1121f')

ggplot(data = pca_data, 
       aes(x = PC1, 
           y = PC2))+
  geom_point(aes(col = Ecotype, 
                 shape = POP_PAIR), 
             size = 2)+
  scale_color_manual(values = plot_pal)


pca_data %>% 
  group_by(POP_PAIR, 
           Ecotype) %>% 
  summarize(n = n())

pca_data %>% 
  filter(POP_PAIR == 'ASHN') %>% 
ggplot(aes(x = PC1, 
           y = PC2))+
  geom_point(aes(col = Ecotype, 
                 shape = POP_PAIR), 
             size = 2)+
  scale_color_manual(values = plot_pal)

pca_data %>% 
  filter(POP_PAIR == 'MYV') %>% 
  ggplot(aes(x = PC1, 
             y = PC2))+
  geom_point(aes(col = Ecotype, 
                 shape = POP_PAIR), 
             size = 2)+
  scale_color_manual(values = plot_pal)

pca_data %>% 
  filter(POP_PAIR == 'SKR') %>% 
  ggplot(aes(x = PC1, 
             y = PC2))+
  geom_point(aes(col = Ecotype, 
                 shape = POP_PAIR), 
             size = 2)+
  scale_color_manual(values = plot_pal)

pca_data %>% 
  filter(POP_PAIR == 'GTS-GAR') %>% 
  ggplot(aes(x = PC1, 
             y = PC2))+
  geom_point(aes(col = Ecotype, 
                 shape = POP_PAIR), 
             size = 2)+
  scale_color_manual(values = plot_pal)


## As sign doesn't matter in a PCA, flip ASHN and GTS axis 1


# Kmeans clustering inversion karyotype -----------------------------------


flip_coords = pca_data %>% 
  filter(POP_PAIR %in% c('ASHN', 
                         'GTS-GAR')) %>% 
  dplyr::select(1:5, 
                Ecotype, 
                POP_PAIR)
  
flip_coords$PC1 = -flip_coords$PC1
flip_coords$PC2 = -flip_coords$PC2

proper_coords = pca_data %>% 
  filter(POP_PAIR %in% c('MYV', 
                         'SKR')) %>% 
  dplyr::select(1:5, 
                Ecotype, 
                POP_PAIR)

clean_pca_data = bind_rows(proper_coords, 
                           flip_coords)




set.seed(1738)

kmean_data = clean_pca_data %>% 
  dplyr::select(PC1, 
                PC2)
km_clust = kmeans(kmean_data, 
                centers = 3, 
                nstart = 25)
kmeans_clusters = km_clust$cluster %>% 
  as.data.frame() %>% 
  as_tibble() %>% 
  rename(kmean_cluster = 1)


cleaner_pca_data = bind_cols(clean_pca_data, 
          kmeans_clusters)



# Final PCA plot - chrxxi inversion ---------------------------------------


chrxxi_inversion_pca = cleaner_pca_data %>% 
  ggplot(aes(x = PC1, 
             y = PC2, 
             group = kmean_cluster))+
  geom_point(aes(col = Ecotype, 
                 shape = POP_PAIR), 
             size = 3)+
  stat_ellipse(type = 'norm')+
  scale_color_manual(values = plot_pal)+
  scale_shape_manual(values = c(15, 16, 17, 18))+
  labs(x = 'Principal component 1 (55.4%)', 
       y = 'Principal component 2 (3.5%)')+
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12), 
        legend.position = 'none')

ggsave(filename = 'ChrXXI_Ecotype_Inversion_plot.tiff', 
       plot = chrxxi_inversion_pca, 
       dpi = 'retina', 
       units = 'cm')




# chrxxi karyotype summary stats ------------------------------------------
cleaner_pca_data %>% 
  group_by(POP_PAIR,
           Ecotype,
           kmean_cluster) %>% 
  summarize(n = n())


cleaner_pca_data %>% 
  group_by(Ecotype,
           kmean_cluster) %>% 
  summarize(n = n())

22/54*100
21/54*100
11/54*100

3/55*100
27/55*100
25/55*100



# heterozygosity per cluster ----------------------------------------------

chrxxi_dart = dartR.base::gl.read.vcf('chrxxi_vcf.vcf')

chrxxi_dart@other$loc.metrics$position = chrxxi_dart@position

filtered_gl_data = dartR.base::gl.filter.locmetric(
  x = chrxxi_dart,
  metric = "position", # The name of the metric you added
  lower = 9963830,
  upper = 11574024,
  keep = "within",
  verbose = 3 # For detailed output
)

ind_names = chrxxi_dart@ind.names %>% 
  as.data.frame() %>% 
  as_tibble() %>% 
  rename(ind.names = 1)

cleaner_pca_data %>% 
  select(Ecotype, 
         POP_PAIR, 
         population, 
         kmean_cluster) %>% 
  bind_cols(., 
            ind_names) %>% 
  dplyr::select(ind.names, 
                # Ecotype, 
                # POP_PAIR, 
                # population, 
                kmean_cluster) %>% 
  rename(id = ind.names, 
         pop = kmean_cluster) %>% 
  write_csv('chrxxi_inversion_metadata.csv')


filtered_gl_data = dartR.base::gl.add.indmetrics(filtered_gl_data, 
                                        ind.metafile = "chrxxi_inversion_metadata.csv")

filtered_gl_data$pop

het_results = dartR.base::gl.report.heterozygosity(filtered_gl_data, 
                                       method = "pop")


het_results$FIS
het_results$Ho
het_results$He

het_data = het_results %>% 
  as.data.frame() %>% 
  as_tibble() %>% 
  mutate(.data = .,
         Karyotype = as.factor(case_when(
           pop == '1' ~ '0/1',
           pop == '2' ~ '1/1',
           pop == '3' ~ '0/0')))

box_plot_pal = c('#003049', 
                 '#7678ed',
                 '#c1121f')

het_data$pop <- factor(het_data$pop,
                       levels = c("3", "1", "2"))
het_data$Karyotype <- factor(het_data$Karyotype,
                       levels = c("0/0", "0/1", "1/1"))


Karyotype_FIS_plot = ggplot(data = het_data, 
       aes(x = Karyotype, 
             y = FIS)) +
  geom_boxplot(aes(col = pop, 
                   fill = pop))+
  scale_fill_manual(values = box_plot_pal)+
  scale_color_manual(values = box_plot_pal)+
  geom_errorbar(aes(ymax = FIS + FISSD, 
                    ymin = FIS - FISSD),
                position = "dodge")+
  theme(panel.grid = element_blank(), 
        axis.title.y = element_text(size = 14), 
        axis.title.x = element_blank(),
        axis.text = element_text(size = 12), 
        legend.position = 'none')

Karyotype_Ho_plot = ggplot(data = het_data, 
       aes(x = Karyotype, 
           y = Ho)) +
  geom_boxplot(aes(col = pop, 
                   fill = pop))+
  scale_fill_manual(values = box_plot_pal)+
  scale_color_manual(values = box_plot_pal)+
  geom_errorbar(aes(ymax = Ho + HoSD, 
                    ymin = Ho - HoSD),
                position = "dodge")+
  theme(panel.grid = element_blank(), 
        axis.title.y = element_text(size = 14), 
        axis.title.x = element_blank(),
        axis.text = element_text(size = 12), 
        legend.position = 'none')


inversion_combo_plots = chrxxi_inversion_pca/(Karyotype_FIS_plot+Karyotype_Ho_plot)

ggsave(filename = 'chrxxi_inversion_plots.tiff', 
       plot = inversion_combo_plots, 
       dpi = 'retina', 
       units = 'cm', 
       width = 30, 
       height = 15)

# PCA per AFVAPER window --------------------------------------------------


# window 1 ----------------------------------------------------------------

win1 = dartR.base::gl.filter.locmetric(
  x = chrxxi_dart,
  metric = "position", # The name of the metric you added
  lower = 9963830,
  upper = 10100975,
  keep = "within",
  verbose = 3 # For detailed output
)

pca_win1 = adegenet::glPca(win1, 
                                nf = 10)

win1_pca_load = pca_win1$loadings %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  as.tibble() %>% 
  rename(SNP = rowname)

win1_pca_score = pca_win1$scores %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  as.tibble() %>% 
  rename(Individual = rowname) 


metadata = read_csv('stickleback_identifiers.csv')
win1_pca_data = bind_cols(metadata, 
                     win1_pca_score) %>% 
  mutate(.data = .,
         Ecotype = as.factor(case_when(
           population == 'ASHNC' ~ 'Cold',
           population == 'ASHNW' ~ 'Warm',
           population == 'CSWY' ~ 'Cold',
           population == 'GTS' ~ 'Warm',
           population == 'MYVC' ~ 'Cold',
           population == 'MYVW' ~ 'Warm',
           population == 'SKRC' ~ 'Cold',
           population == 'SKRW' ~ 'Warm')))%>% 
  mutate(.data = .,
         POP_PAIR = as.factor(case_when(
           population == 'ASHNC' ~ 'ASHN',
           population == 'ASHNW' ~ 'ASHN',
           population == 'CSWY' ~ 'GTS-GAR',
           population == 'GTS' ~ 'GTS-GAR',
           population == 'MYVC' ~ 'MYV',
           population == 'MYVW' ~ 'MYV',
           population == 'SKRC' ~ 'SKR',
           population == 'SKRW' ~ 'SKR')))

plot_pal = c('#003049', 
             '#c1121f')

ggplot(data = win1_pca_data, 
       aes(x = PC1, 
           y = PC2))+
  geom_point(aes(col = Ecotype, 
                 shape = POP_PAIR), 
             size = 2)+
  scale_color_manual(values = plot_pal)+
  labs(title = 'Window 1: 9963830-10100975')



# window 2 ----------------------------------------------------------------

win2 = dartR.base::gl.filter.locmetric(
  x = chrxxi_dart,
  metric = "position", # The name of the metric you added
  lower = 10100983,
  upper = 10264227,
  keep = "within",
  verbose = 3 # For detailed output
)

pca_win2 = adegenet::glPca(win2, 
                           nf = 10)

win2_pca_load = pca_win2$loadings %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  as.tibble() %>% 
  rename(SNP = rowname)

win2_pca_score = pca_win2$scores %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  as.tibble() %>% 
  rename(Individual = rowname) 


metadata = read_csv('stickleback_identifiers.csv')
win2_pca_data = bind_cols(metadata, 
                          win2_pca_score) %>% 
  mutate(.data = .,
         Ecotype = as.factor(case_when(
           population == 'ASHNC' ~ 'Cold',
           population == 'ASHNW' ~ 'Warm',
           population == 'CSWY' ~ 'Cold',
           population == 'GTS' ~ 'Warm',
           population == 'MYVC' ~ 'Cold',
           population == 'MYVW' ~ 'Warm',
           population == 'SKRC' ~ 'Cold',
           population == 'SKRW' ~ 'Warm')))%>% 
  mutate(.data = .,
         POP_PAIR = as.factor(case_when(
           population == 'ASHNC' ~ 'ASHN',
           population == 'ASHNW' ~ 'ASHN',
           population == 'CSWY' ~ 'GTS-GAR',
           population == 'GTS' ~ 'GTS-GAR',
           population == 'MYVC' ~ 'MYV',
           population == 'MYVW' ~ 'MYV',
           population == 'SKRC' ~ 'SKR',
           population == 'SKRW' ~ 'SKR')))

plot_pal = c('#003049', 
             '#c1121f')

ggplot(data = win2_pca_data, 
       aes(x = PC1, 
           y = PC2))+
  geom_point(aes(col = Ecotype, 
                 shape = POP_PAIR), 
             size = 2)+
  scale_color_manual(values = plot_pal)+
  labs(title = 'Window 2: 10100983-10264227')



# window 3 ----------------------------------------------------------------

win3 = dartR.base::gl.filter.locmetric(
  x = chrxxi_dart,
  metric = "position", # The name of the metric you added
  lower = 10266044,
  upper = 10403157,
  keep = "within",
  verbose = 3 # For detailed output
)

pca_win3 = adegenet::glPca(win3, 
                           nf = 10)

win3_pca_load = pca_win3$loadings %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  as.tibble() %>% 
  rename(SNP = rowname)

win3_pca_score = pca_win3$scores %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  as.tibble() %>% 
  rename(Individual = rowname) 


metadata = read_csv('stickleback_identifiers.csv')
win3_pca_data = bind_cols(metadata, 
                          win3_pca_score) %>% 
  mutate(.data = .,
         Ecotype = as.factor(case_when(
           population == 'ASHNC' ~ 'Cold',
           population == 'ASHNW' ~ 'Warm',
           population == 'CSWY' ~ 'Cold',
           population == 'GTS' ~ 'Warm',
           population == 'MYVC' ~ 'Cold',
           population == 'MYVW' ~ 'Warm',
           population == 'SKRC' ~ 'Cold',
           population == 'SKRW' ~ 'Warm')))%>% 
  mutate(.data = .,
         POP_PAIR = as.factor(case_when(
           population == 'ASHNC' ~ 'ASHN',
           population == 'ASHNW' ~ 'ASHN',
           population == 'CSWY' ~ 'GTS-GAR',
           population == 'GTS' ~ 'GTS-GAR',
           population == 'MYVC' ~ 'MYV',
           population == 'MYVW' ~ 'MYV',
           population == 'SKRC' ~ 'SKR',
           population == 'SKRW' ~ 'SKR')))

plot_pal = c('#003049', 
             '#c1121f')

ggplot(data = win3_pca_data, 
       aes(x = PC1, 
           y = PC2))+
  geom_point(aes(col = Ecotype, 
                 shape = POP_PAIR), 
             size = 2)+
  scale_color_manual(values = plot_pal)+
  labs(title = 'window 3: 10266044-10403157')


