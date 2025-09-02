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


