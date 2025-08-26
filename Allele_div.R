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

vcf = read.vcfR('stickleback_filtered_vcf.vcf')

head(vcf)

chrxxi_data = vcf[getCHROM(vcf) == 'chr_XXI']

head(chrxxi_data)

chrxxi_data@gt


chrxxi_data@gt %>% 
  as.data.frame() %>% 
  as_tibble()
