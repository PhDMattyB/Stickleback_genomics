##############################
## Temperature associated outliers
##
## Matt Brachmann (PhDMattyB)
##
## 01.09.2025
##
##############################

theme_set(theme_bw())

# setwd('C:/Stickleback_Genomic/vcf_filter/')
setwd('~/Parsons_Postdoc/Stickleback_Genomic/vcf_filter/')




##Common FST outliers
WC_Fst_clean_outs = read_csv('WC_Fst_clean.csv') %>% 
  stickle_CHR_reorder() %>% 
  dist_cal()%>% 
  filter(value == 'Outlier')


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



FST_outs_LFMM %>% 
  filter(CHR == 'chr_XXI') %>% View()
