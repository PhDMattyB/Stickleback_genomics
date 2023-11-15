##############################
## selective sweeps iHS XP-EHH
##
## Matt Brachmann (PhDMattyB)
##
## 10.10.2023
##
##############################


setwd('~/Parsons_Postdoc/Stickleback_Genomic/vcf_filter/')

library(tidyverse)
# install.packages('rehh')
library(rehh)

theme_set(theme_bw())

#
# MYVATN signatures of selection ------------------------------------------


myvc = data2haplohh(hap_file = 'MYVC_genotypes.vcf', 
                    polarize_vcf = F, 
                    min_maf = 0.05, 
                    chr.name = 'chr_XXI', 
                    allele_coding = 'map')

myvw = data2haplohh(hap_file = 'MYVW_genotypes.vcf', 
                    polarize_vcf = F, 
                    min_maf = 0.05, 
                    chr.name = 'chr_XXI', 
                    allele_coding = 'map')


# iHS haplotype scan

myvc_scan = scan_hh(myvc, polarized = F)
myvw_scan = scan_hh(myvw, polarized = F)

## perform ihs selective scan
myvw_ihs = ihh2ihs(myvw_scan, 
                   freqbin = 1)
myvc_ihs = ihh2ihs(myvc_scan, 
                   freqbin = 1)

## plot the ihs statistic
ggplot(myvw_ihs$ihs, 
       aes(POSITION, 
           IHS))+
  geom_point()

ggplot(myvc_ihs$ihs, 
       aes(POSITION, 
           IHS))+
  geom_point()

## plot the pvalues
ggplot(myvw_ihs$ihs, 
       aes(POSITION, 
           LOGPVALUE))+
  geom_point()

ggplot(myvc_ihs$ihs, 
       aes(POSITION, 
           LOGPVALUE))+
  geom_point()



# xp-ehh analysis (cross population) 
myvw_myvc = ies2xpehh(myvw_scan, 
                      myvc_scan, 
                      popname1 = 'MYVW', 
                      popname2 = 'MYVC', 
                      include_freq = T)

# plot

ggplot(myvw_myvc, 
       aes(POSITION, 
           XPEHH_MYVW_MYVC))+
  geom_point()

ggplot(myvw_myvc, 
       aes(POSITION, 
           LOGPVALUE)) + geom_point()



# Haplotype structure around selection target 
hit = myvw_myvc %>% 
  arrange(desc(LOGPVALUE)) %>% 
  top_n(1)

# get SNP position
x = hit$POSITION

marker_id_warm = which(myvw@positions == x)
marker_id_cold = which(myvc@positions == x)

myvw_furcation = calc_furcation(myvw, 
                                mrk = marker_id_warm)

myvc_furcation = calc_furcation(myvc, 
                                mrk = marker_id_cold )


plot(myvw_furcation, 
     xlim = c(58695, 17420697))
plot(myvc_furcation, 
     xlim = c(58695, 17420697))

myvw_haplen = calc_haplen(myvw_furcation)
myvc_haplen = calc_haplen(myvc_furcation)

plot(myvw_haplen)
plot(myvc_haplen)


# ASHN signatures of selection ------------------------------------------


ashnc = data2haplohh(hap_file = 'ashnC_genotypes.vcf', 
                    polarize_vcf = F, 
                    min_maf = 0.05, 
                    chr.name = 'chr_XXI', 
                    allele_coding = 'map')

ashnw = data2haplohh(hap_file = 'ashnW_genotypes.vcf', 
                    polarize_vcf = F, 
                    min_maf = 0.05, 
                    chr.name = 'chr_XXI', 
                    allele_coding = 'map')


# iHS haplotype scan

ashnc_scan = scan_hh(ashnc, polarized = F)
ashnw_scan = scan_hh(ashnw, polarized = F)

## perform ihs selective scan
ashnw_ihs = ihh2ihs(ashnw_scan, 
                   freqbin = 1)
ashnc_ihs = ihh2ihs(ashnc_scan, 
                   freqbin = 1)

## plot the ihs statistic
ggplot(ashnw_ihs$ihs, 
       aes(POSITION, 
           IHS))+
  geom_point()

ggplot(ashnc_ihs$ihs, 
       aes(POSITION, 
           IHS))+
  geom_point()

## plot the pvalues
ggplot(ashnw_ihs$ihs, 
       aes(POSITION, 
           LOGPVALUE))+
  geom_point()

ggplot(ashnc_ihs$ihs, 
       aes(POSITION, 
           LOGPVALUE))+
  geom_point()



# xp-ehh analysis (cross population) 
ashnw_ashnc = ies2xpehh(ashnw_scan, 
                      ashnc_scan, 
                      popname1 = 'ashnW', 
                      popname2 = 'ashnC', 
                      include_freq = T)

# plot

ggplot(ashnw_ashnc, 
       aes(POSITION, 
           XPEHH_ashnW_ashnC))+
  geom_point()

ggplot(ashnw_ashnc, 
       aes(POSITION, 
           LOGPVALUE)) + geom_point()



# Haplotype structure around selection target 
# hit = ashnw_ashnc %>% 
#   arrange(desc(LOGPVALUE)) %>% 
#   top_n(1)

# get SNP position
# x = hit$POSITION
# 
# marker_id_warm = which(ashnw@positions == x)
# marker_id_cold = which(ashnc@positions == x)
# 
# ashnw_furcation = calc_furcation(ashnw, 
#                                 mrk = marker_id_warm)
# 
# ashnc_furcation = calc_furcation(ashnc, 
#                                 mrk = marker_id_cold )
# 
# 
# plot(ashnw_furcation, 
#      xlim = c(58695, 17420697))
# plot(ashnc_furcation, 
#      xlim = c(58695, 17420697))
# 
# ashnw_haplen = calc_haplen(ashnw_furcation)
# ashnc_haplen = calc_haplen(ashnc_furcation)
# 
# plot(ashnw_haplen)
# plot(ashnc_haplen)

# SKR signatures of selection ------------------------------------------


skrc = data2haplohh(hap_file = 'skrC_genotypes.vcf', 
                    polarize_vcf = F, 
                    min_maf = 0.05, 
                    chr.name = 'chr_XXI', 
                    allele_coding = 'map')

skrw = data2haplohh(hap_file = 'skrW_genotypes.vcf', 
                    polarize_vcf = F, 
                    min_maf = 0.05, 
                    chr.name = 'chr_XXI', 
                    allele_coding = 'map')


# iHS haplotype scan

skrc_scan = scan_hh(skrc, polarized = F)
skrw_scan = scan_hh(skrw, polarized = F)

## perform ihs selective scan
skrw_ihs = ihh2ihs(skrw_scan, 
                   freqbin = 1)
skrc_ihs = ihh2ihs(skrc_scan, 
                   freqbin = 1)

## plot the ihs statistic
ggplot(skrw_ihs$ihs, 
       aes(POSITION, 
           IHS))+
  geom_point()

ggplot(skrc_ihs$ihs, 
       aes(POSITION, 
           IHS))+
  geom_point()

## plot the pvalues
ggplot(skrw_ihs$ihs, 
       aes(POSITION, 
           LOGPVALUE))+
  geom_point()

ggplot(skrc_ihs$ihs, 
       aes(POSITION, 
           LOGPVALUE))+
  geom_point()



# xp-ehh analysis (cross population) 
skrw_skrc = ies2xpehh(skrw_scan, 
                      skrc_scan, 
                      popname1 = 'skrW', 
                      popname2 = 'skrC', 
                      include_freq = T)

# plot

ggplot(skrw_skrc, 
       aes(POSITION, 
           XPEHH_skrW_skrC))+
  geom_point()

ggplot(skrw_skrc, 
       aes(POSITION, 
           LOGPVALUE)) + geom_point()



# Haplotype structure around selection target 
# hit = skrw_skrc %>% 
#   arrange(desc(LOGPVALUE)) %>% 
#   top_n(1)
# 
# # get SNP position
# x = hit$POSITION
# 
# marker_id_warm = which(skrw@positions == x)
# marker_id_cold = which(skrc@positions == x)
# 
# skrw_furcation = calc_furcation(skrw, 
#                                 mrk = marker_id_warm)
# 
# skrc_furcation = calc_furcation(skrc, 
#                                 mrk = marker_id_cold )
# 
# 
# plot(skrw_furcation, 
#      xlim = c(58695, 17420697))
# plot(skrc_furcation, 
#      xlim = c(58695, 17420697))
# 
# skrw_haplen = calc_haplen(skrw_furcation)
# skrc_haplen = calc_haplen(skrc_furcation)
# 
# plot(skrw_haplen)
# plot(skrc_haplen)

# GTS signatures of selection ------------------------------------------


gts = data2haplohh(hap_file = 'GTS_genotypes.vcf', 
                    polarize_vcf = F, 
                    min_maf = 0.05, 
                    chr.name = 'chr_XXI', 
                    allele_coding = 'map')


# iHS haplotype scan

gts_scan = scan_hh(gts, polarized = F)

## perform ihs selective scan
gts_ihs = ihh2ihs(gts_scan, 
                   freqbin = 1)

## plot the ihs statistic
ggplot(gts_ihs$ihs, 
       aes(POSITION, 
           IHS))+
  geom_point()


## plot the pvalues
ggplot(gts_ihs$ihs, 
       aes(POSITION, 
           LOGPVALUE))+
  geom_point()


# CSWY signatures of selection ------------------------------------------

cswy = data2haplohh(hap_file = 'CSWY_genotypes.vcf', 
                    polarize_vcf = F, 
                    min_maf = 0.05, 
                    chr.name = 'chr_XXI', 
                    allele_coding = 'map')

# iHS haplotype scan

cswy_scan = scan_hh(cswy, polarized = F)

## perform ihs selective scan
cswy_ihs = ihh2ihs(cswy_scan, 
                   freqbin = 1)

## plot the ihs statistic
ggplot(cswy_ihs$ihs, 
       aes(POSITION, 
           IHS))+
  geom_point()


## plot the pvalues
ggplot(cswy_ihs$ihs, 
       aes(POSITION, 
           LOGPVALUE))+
  geom_point()


# Warm vs Cold sweep detect (chr 21) --------------------------------------

warm = data2haplohh(hap_file = 'Warm_morphs_genotypes.vcf', 
                    polarize_vcf = F, 
                    min_maf = 0.05, 
                    chr.name = 'chr_XXI', 
                    allele_coding = 'map')

cold = data2haplohh(hap_file = 'Cold_morphs_genotypes.vcf', 
                    polarize_vcf = F, 
                    min_maf = 0.05, 
                    chr.name = 'chr_XXI', 
                    allele_coding = 'map')


# iHS haplotype scan 


cold_scan = scan_hh(cold, polarized = F)
warm_scan = scan_hh(warm, polarized = F)

## perform ihs selective scan
warm_ihs = ihh2ihs(warm_scan, 
                   freqbin = 1)
cold_ihs = ihh2ihs(cold_scan, 
                   freqbin = 1)

## plot the ihs statistic
ggplot(warm_ihs$ihs, 
       aes(POSITION, 
           IHS))+
  geom_point()

ggplot(cold_ihs$ihs, 
       aes(POSITION, 
           IHS))+
  geom_point()

## plot the pvalues
ggplot(warm_ihs$ihs, 
       aes(POSITION, 
           LOGPVALUE))+
  geom_point()

ggplot(cold_ihs$ihs, 
       aes(POSITION, 
           LOGPVALUE))+
  geom_point()



# xp-ehh analysis (cross population) 

warm_cold = ies2xpehh(warm_scan, 
                      cold_scan, 
                      popname1 = 'warm', 
                      popname2 = 'cold', 
                      include_freq = T)

# plot

ggplot(warm_cold, 
       aes(POSITION, 
           XPEHH_warm_cold))+
  geom_point()

ggplot(warm_cold, 
       aes(POSITION, 
           LOGPVALUE)) + geom_point()



# Haplotype structure around selection target 
# # find the highest hit
# hit = warm_cold %>% 
#   arrange(desc(LOGPVALUE)) %>% 
#   top_n(1)
# 
# # get SNP position
# x = hit$POSITION
# 
# marker_id_warm = which(warm@positions == x)
# marker_id_cold = which(cold@positions == x)
# 
# warm_furcation = calc_furcation(warm, 
#                                 mrk = marker_id_warm)
# 
# cold_furcation = calc_furcation(cold, 
#                                 mrk = marker_id_cold )
# 
# 
# plot(warm_furcation, 
#      xlim = c(58695, 17420697))
# plot(cold_furcation, 
#      xlim = c(58694, 17420697))
# 
# warm_haplen = calc_haplen(warm_furcation)
# cold_haplen = calc_haplen(cold_furcation)
# 
# plot(warm_haplen)
# plot(cold_haplen)
# 

# Warm vs Cold sweep detect (chr 21) NO GTS or CSWY --------------------------------------

warm2 = data2haplohh(hap_file = 'Warm_morphs_NOGTS_genotypes.vcf', 
                    polarize_vcf = F, 
                    min_maf = 0.05, 
                    chr.name = 'chr_XXI', 
                    allele_coding = 'map')


cold2 = data2haplohh(hap_file = 'Cold_morphs_NOCSWY_genotypes.vcf', 
                    polarize_vcf = F, 
                    min_maf = 0.05, 
                    chr.name = 'chr_XXI', 
                    allele_coding = 'map')


# iHS haplotype scan 


cold2_scan = scan_hh(cold2, polarized = F)
warm2_scan = scan_hh(warm2, polarized = F)

## perform ihs selective scan
warm2_ihs = ihh2ihs(warm2_scan, 
                   freqbin = 1)
cold2_ihs = ihh2ihs(cold2_scan, 
                   freqbin = 1)

# ## plot the ihs statistic
# ggplot(warm2_ihs$ihs, 
#        aes(POSITION, 
#            IHS))+
#   geom_point()
# 
# ggplot(cold2_ihs$ihs, 
#        aes(POSITION, 
#            IHS))+
#   geom_point()
# 
# ## plot the pvalues
# ggplot(warm2_ihs$ihs, 
#        aes(POSITION, 
#            LOGPVALUE))+
#   geom_point()
# 
# ggplot(cold2_ihs$ihs, 
#        aes(POSITION, 
#            LOGPVALUE))+
#   geom_point()



# xp-ehh analysis (cross population) 

warm2_cold2 = ies2xpehh(warm2_scan, 
                      cold2_scan, 
                      popname1 = 'warm2', 
                      popname2 = 'cold2', 
                      include_freq = T,
                      p.adjust.method = 'bonferroni')

calc_candidate_regions(warm2_cold2,
                       threshold = 2.313087,
                       pval = TRUE,
                       window_size = 1E6,
                       overlap = 1E5,
                       min_n_extr_mrk = 2)


# plot

warm_cold_xpehh = ggplot(warm2_cold2, 
       aes(POSITION, 
           XPEHH_warm2_cold2))+
  geom_point()+
  geom_hline(yintercept = 0, 
             col = '#d00000', 
             linewidth = 1)+
  # ylim(-10, 10)+ 
  labs(x = 'Position (bp)', 
       y = 'XP-EHH statistic')+
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12))

ggsave('~/Parsons_Postdoc/Stickleback_Genomic/Figures/warm_cold_xpehh_stat.tiff', 
       plot = warm_cold_xpehh, 
       dpi = 'retina', 
       units = 'cm',
       width = 15, 
       height = 10)


warm_cold_xpehh_pval = ggplot(warm2_cold2, 
       aes(POSITION, 
           LOGPVALUE))+
  geom_point()+
  geom_hline(yintercept = 2.313087, 
             col = '#d00000', 
             linewidth = 1)+
  # ylim(-10, 10)+ 
  labs(x = 'Position (bp)', 
       y = 'Log p-value')+
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12))

ggsave('~/Parsons_Postdoc/Stickleback_Genomic/Figures/warm_cold_xpehh_pvalues.tiff', 
       plot = warm_cold_xpehh_pval, 
       dpi = 'retina', 
       units = 'cm',
       width = 15, 
       height = 10)



unphased_sweep_outliers = warm2_cold2 %>% 
  as_tibble() %>% 
  dplyr::select(CHR, 
                POSITION, 
                XPEHH_warm2_cold2, 
                LOGPVALUE) %>% 
  mutate(pvalue = 10**-LOGPVALUE) %>% 
  na.omit() %>% 
  # filter(LOGPVALUE >= 1.39)
  filter(pvalue <= 0.05)

View(unphased_sweep_outliers)


# Haplotype structure around selection target 
# find the highest hit
# hit = warm2_cold2 %>% 
#   arrange(desc(LOGPVALUE)) %>% 
#   top_n(10)
# 
# 
# ## need to center this around the peak of the inversion site
# ## Right now it's picking out the tip of a separate selective sweep
# # get SNP position
# x = hit$POSITION
# 
# marker_id_warm2 = which(warm2@positions == x)
# marker_id_cold2 = which(cold2@positions == x)
# 
# warm2_furcation = calc_furcation(warm2, 
#                                 mrk = marker_id_warm2)
# 
# cold2_furcation = calc_furcation(cold2, 
#                                 mrk = marker_id_cold2 )
# 
# 
# plot(warm2_furcation, 
#      xlim = c(58695, 17420697))
# plot(cold2_furcation, 
#      xlim = c(58694, 17420697))
# 
# warm2_haplen = calc_haplen(warm2_furcation)
# cold2_haplen = calc_haplen(cold2_furcation)
# 
# plot(warm2_haplen)
# plot(cold2_haplen)





# Phased analyses ---------------------------------------------------------

# Warm vs Cold sweep detect (chr 21) NO GTS or CSWY --------------------------------------

warm_phased = data2haplohh(hap_file = 'Warm_morphs_NOGTS_phased_chr21', 
                     polarize_vcf = F, 
                     min_maf = 0.05, 
                     chr.name = 'chr_XXI', 
                     allele_coding = 'map')


cold_phased = data2haplohh(hap_file = 'Cold_morphs_NOCSWY_phased_chr21', 
                     polarize_vcf = F, 
                     min_maf = 0.05, 
                     chr.name = 'chr_XXI', 
                     allele_coding = 'map')


# iHS haplotype scan 


cold_phased_scan = scan_hh(cold_phased, polarized = F)
warm_phased_scan = scan_hh(warm_phased, polarized = F)

## perform ihs selective scan
warm_phased_ihs = ihh2ihs(warm_phased_scan, 
                    freqbin = 1)
cold_phased_ihs = ihh2ihs(cold_phased_scan, 
                    freqbin = 1)

# ## plot the ihs statistic
# ggplot(warm_phased_ihs$ihs, 
#        aes(POSITION, 
#            IHS))+
#   geom_point()
# 
# ggplot(cold_phased_ihs$ihs, 
#        aes(POSITION, 
#            IHS))+
#   geom_point()
# 
# ## plot the pvalues
# ggplot(warm_phased_ihs$ihs, 
#        aes(POSITION, 
#            LOGPVALUE))+
#   geom_point()
# 
# ggplot(cold_phased_ihs$ihs, 
#        aes(POSITION, 
#            LOGPVALUE))+
#   geom_point()



# xp-ehh analysis (cross population) 

warm_phased_cold_phased = ies2xpehh(warm_phased_scan, 
                        cold_phased_scan, 
                        popname1 = 'warm_phased', 
                        popname2 = 'cold_phased', 
                        include_freq = T, 
                        p.adjust.method = 'bonferroni')

# plot

warm_cold_xpehh_phased = ggplot(warm_phased_cold_phased, 
                         aes(POSITION, 
                             XPEHH_warm_phased_cold_phased))+
  geom_point()+
  geom_hline(yintercept = 0, 
             col = '#2a9d8f', 
             linewidth = 1)+
  # ylim(-10, 10)+ 
  labs(x = 'Position (bp)', 
       y = 'XP-EHH statistic')+
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12))

ggsave('~/Parsons_Postdoc/Stickleback_Genomic/Figures/warm_cold_xpehh_stat_PHASED.tiff', 
       plot = warm_cold_xpehh_phased, 
       dpi = 'retina', 
       units = 'cm',
       width = 20.0, 
       height = 13)


warm_cold_xpehh_pvalues = ggplot(warm_phased_cold_phased, 
                                 aes(POSITION, 
                                     LOGPVALUE)) + 
  geom_point()+
  geom_hline(yintercept = 5, 
             col = '#2a9d8f', 
             linewidth = 1)+
  labs(x = 'Position (bp)', 
       y = 'Log p-value')+
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12))

ggsave('~/Parsons_Postdoc/Stickleback_Genomic/Figures/warm_cold_xpehh_pvalues_PHASED.tiff', 
       plot = warm_cold_xpehh_pvalues, 
       dpi = 'retina', 
       units = 'cm',
       width = 20.0, 
       height = 13)

phased_sweep_outliers = warm_phased_cold_phased %>% 
  as_tibble() %>% 
  dplyr::select(CHR, 
                POSITION, 
                XPEHH_warm_phased_cold_phased, 
                LOGPVALUE) %>% 
  mutate(pvalue = 10**-LOGPVALUE) %>% 
  na.omit() %>% 
  filter(pvalue <= 0.05)

# Haplotype structure around selection target 
# find the highest hit
# hit = warm_phased_cold_phased %>% 
#   arrange(desc(LOGPVALUE)) %>% 
#   top_n(10)
# 
# 
# ## need to center this around the peak of the inversion site
# ## Right now it's picking out the tip of a separate selective sweep
# # get SNP position
# x = hit$POSITION
# 
# marker_id_warm_phased = which(warm_phased@positions == x)
# marker_id_cold_phased = which(cold_phased@positions == x)
# 
# warm_phased_furcation = calc_furcation(warm_phased, 
#                                 mrk = marker_id_warm_phased)
# 
# cold_phased_furcation = calc_furcation(cold_phased, 
#                                 mrk = marker_id_cold_phased )
# 
# 
# plot(warm_phased_furcation, 
#      xlim = c(58695, 17420697))
# plot(cold_phased_furcation, 
#      xlim = c(58694, 17420697))
# 
# warm_phased_haplen = calc_haplen(warm_phased_furcation)
# cold_phased_haplen = calc_haplen(cold_phased_furcation)
# 
# plot(warm_phased_haplen)
# plot(cold_phased_haplen)

# Warm vs Cold sweep detect (chr 21) NO GTS or CSWY window 0.5 --------------------------------------

warm_0.5 = data2haplohh(hap_file = 'Warm_morphs_NOGTS_phased_chr21_window0.5', 
                     polarize_vcf = F, 
                     min_maf = 0.05, 
                     chr.name = 'chr_XXI', 
                     allele_coding = 'map')


cold_0.5 = data2haplohh(hap_file = 'Cold_morphs_NOCSWY_phased_chr21_window0.5', 
                     polarize_vcf = F, 
                     min_maf = 0.05, 
                     chr.name = 'chr_XXI', 
                     allele_coding = 'map')


# iHS haplotype scan 


cold_0.5_scan = scan_hh(cold_0.5, polarized = F)
warm_0.5_scan = scan_hh(warm_0.5, polarized = F)

## perform ihs selective scan
warm_0.5_ihs = ihh2ihs(warm_0.5_scan, 
                    freqbin = 1)
cold_0.5_ihs = ihh2ihs(cold_0.5_scan, 
                    freqbin = 1)

# ## plot the ihs statistic
# ggplot(warm_0.5_ihs$ihs, 
#        aes(POSITION, 
#            IHS))+
#   geom_point()
# 
# ggplot(cold_0.5_ihs$ihs, 
#        aes(POSITION, 
#            IHS))+
#   geom_point()
# 
# ## plot the pvalues
# ggplot(warm_0.5_ihs$ihs, 
#        aes(POSITION, 
#            LOGPVALUE))+
#   geom_point()
# 
# ggplot(cold_0.5_ihs$ihs, 
#        aes(POSITION, 
#            LOGPVALUE))+
#   geom_point()



# xp-ehh analysis (cross population) 

warm_0.5_cold_0.5 = ies2xpehh(warm_0.5_scan, 
                        cold_0.5_scan, 
                        popname1 = 'warm_0.5', 
                        popname2 = 'cold_0.5', 
                        include_freq = T)

# plot

warm_cold_xpehh_phased_win0.5 = ggplot(warm_0.5_cold_0.5, 
                                aes(POSITION, 
                                    XPEHH_warm_0.5_cold_0.5))+
  geom_point()+
  geom_hline(yintercept = 0, 
             col = '#2a9d8f', 
             linewidth = 1)+
  # ylim(-10, 10)+ 
  labs(x = 'Position (bp)', 
       y = 'XP-EHH statistic')+
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12))

ggsave('~/Parsons_Postdoc/Stickleback_Genomic/Figures/warm_cold_xpehh_stat_PHASED.tiff', 
       plot = warm_cold_xpehh_phased, 
       dpi = 'retina', 
       units = 'cm',
       width = 20.0, 
       height = 13)


warm_cold_xpehh_pvalues = ggplot(warm_0.5_cold_0.5, 
                                 aes(POSITION, 
                                     LOGPVALUE)) + 
  geom_point()+
  geom_hline(yintercept = 5, 
             col = '#2a9d8f', 
             linewidth = 1)+
  labs(x = 'Position (bp)', 
       y = 'Log p-value')+
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12))

ggsave('~/Parsons_Postdoc/Stickleback_Genomic/Figures/warm_cold_xpehh_pvalues_PHASED.tiff', 
       plot = warm_cold_xpehh_pvalues, 
       dpi = 'retina', 
       units = 'cm',
       width = 20.0, 
       height = 13)

# Haplotype structure around selection target 
# find the highest hit
# hit = warm_0.5_cold_0.5 %>% 
#   arrange(desc(LOGPVALUE)) %>% 
#   top_n(10)
# 
# 
# ## need to center this around the peak of the inversion site
# ## Right now it's picking out the tip of a separate selective sweep
# # get SNP position
# x = hit$POSITION
# 
# marker_id_warm_0.5 = which(warm_0.5@positions == x)
# marker_id_cold_0.5 = which(cold_0.5@positions == x)
# 
# warm_0.5_furcation = calc_furcation(warm_0.5, 
#                                 mrk = marker_id_warm_0.5)
# 
# cold_0.5_furcation = calc_furcation(cold_0.5, 
#                                 mrk = marker_id_cold_0.5 )
# 
# 
# plot(warm_0.5_furcation, 
#      xlim = c(58695, 17420697))
# plot(cold_0.5_furcation, 
#      xlim = c(58694, 17420697))
# 
# warm_0.5_haplen = calc_haplen(warm_0.5_furcation)
# cold_0.5_haplen = calc_haplen(cold_0.5_furcation)
# 
# plot(warm_0.5_haplen)
# plot(cold_0.5_haplen)

# Warm vs Cold sweep detect (chr 21) NO GTS or CSWY window 0.25 --------------------------------------

warm_0.25 = data2haplohh(hap_file = 'Warm_morphs_NOGTS_phased_chr21_window0.25', 
                        polarize_vcf = F, 
                        min_maf = 0.05, 
                        chr.name = 'chr_XXI', 
                        allele_coding = 'map')


cold_0.25 = data2haplohh(hap_file = 'Cold_morphs_NOCSWY_phased_chr21_window0.25', 
                        polarize_vcf = F, 
                        min_maf = 0.05, 
                        chr.name = 'chr_XXI', 
                        allele_coding = 'map')


# iHS haplotype scan 


cold_0.25_scan = scan_hh(cold_0.25, polarized = F)
warm_0.25_scan = scan_hh(warm_0.25, polarized = F)

## perform ihs selective scan
warm_0.25_ihs = ihh2ihs(warm_0.25_scan, 
                       freqbin = 1)
cold_0.25_ihs = ihh2ihs(cold_0.25_scan, 
                       freqbin = 1)

# ## plot the ihs statistic
# ggplot(warm_0.25_ihs$ihs, 
#        aes(POSITION, 
#            IHS))+
#   geom_point()
# 
# ggplot(cold_0.25_ihs$ihs, 
#        aes(POSITION, 
#            IHS))+
#   geom_point()
# 
# ## plot the pvalues
# ggplot(warm_0.25_ihs$ihs, 
#        aes(POSITION, 
#            LOGPVALUE))+
#   geom_point()
# 
# ggplot(cold_0.25_ihs$ihs, 
#        aes(POSITION, 
#            LOGPVALUE))+
#   geom_point()



# xp-ehh analysis (cross population) 

warm_0.25_cold_0.25 = ies2xpehh(warm_0.25_scan, 
                              cold_0.25_scan, 
                              popname1 = 'warm_0.25', 
                              popname2 = 'cold_0.25', 
                              include_freq = T)

# plot

warm_cold_xpehh_phased_win0.25 = ggplot(warm_0.25_cold_0.25, 
                                       aes(POSITION, 
                                           XPEHH_warm_0.25_cold_0.25))+
  geom_point()+
  geom_hline(yintercept = 0, 
             col = '#2a9d8f', 
             linewidth = 1)+
  # ylim(-10, 10)+ 
  labs(x = 'Position (bp)', 
       y = 'XP-EHH statistic')+
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12))

ggsave('~/Parsons_Postdoc/Stickleback_Genomic/Figures/warm_cold_xpehh_stat_PHASED.tiff', 
       plot = warm_cold_xpehh_phased, 
       dpi = 'retina', 
       units = 'cm',
       width = 20.0, 
       height = 13)


warm_cold_xpehh_pvalues = ggplot(warm_0.25_cold_0.25, 
                                 aes(POSITION, 
                                     LOGPVALUE)) + 
  geom_point()+
  geom_hline(yintercept = 5, 
             col = '#2a9d8f', 
             linewidth = 1)+
  labs(x = 'Position (bp)', 
       y = 'Log p-value')+
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12))

ggsave('~/Parsons_Postdoc/Stickleback_Genomic/Figures/warm_cold_xpehh_pvalues_PHASED.tiff', 
       plot = warm_cold_xpehh_pvalues, 
       dpi = 'retina', 
       units = 'cm',
       width = 20.0, 
       height = 13)

# Haplotype structure around selection target 
# find the highest hit
# hit = warm_0.25_cold_0.25 %>% 
#   arrange(desc(LOGPVALUE)) %>% 
#   top_n(10)
# 
# 
# ## need to center this around the peak of the inversion site
# ## Right now it's picking out the tip of a separate selective sweep
# # get SNP position
# x = hit$POSITION
# 
# marker_id_warm_0.25 = which(warm_0.25@positions == x)
# marker_id_cold_0.25 = which(cold_0.25@positions == x)
# 
# warm_0.25_furcation = calc_furcation(warm_0.25, 
#                                 mrk = marker_id_warm_0.25)
# 
# cold_0.25_furcation = calc_furcation(cold_0.25, 
#                                 mrk = marker_id_cold_0.25 )
# 
# 
# plot(warm_0.25_furcation, 
#      xlim = c(58695, 17420697))
# plot(cold_0.25_furcation, 
#      xlim = c(58694, 17420697))
# 
# warm_0.25_haplen = calc_haplen(warm_0.25_furcation)
# cold_0.25_haplen = calc_haplen(cold_0.25_furcation)
# 
# plot(warm_0.25_haplen)
# plot(cold_0.25_haplen)




