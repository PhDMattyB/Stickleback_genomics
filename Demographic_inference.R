##############################
## demographic inference
##
## Matt Brachmann (PhDMattyB)
##
## 23.06.2023
##
##############################

# install.packages('pegas')
# install.packages('adegenet')
library(adegenet)
library(pegas)
library(tidyverse)

setwd('~/Parsons_Postdoc/Stickleback_Genomic/vcf_filter/')

warm_stickle = read.vcf('stickleback_warm_pops.vcf', 
                        from = 1, 
                        to = 173444)

warm_SFS = site.spectrum(warm_stickle, 
              folded = T)

cold_stickle = read.vcf('stickleback_cold_pops.vcf', 
                        from = 1, 
                        to = 173444)
cold_sfs = site.spectrum(cold_stickle, 
                         folded = T)
