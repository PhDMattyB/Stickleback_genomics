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

myvc = data2haplohh(hap_file = 'MYVC_genotypes.vcf', 
                    polarize_vcf = F, 
                    min_maf = 0.05, 
                    chr.name = 'chr_XXI')

myvw = data2haplohh(hap_file = 'MYVW_genotypes.vcf', 
                    polarize_vcf = F, 
                    min_maf = 0.05, 
                    chr.name = 'chr_XXI')


# iHS haplotype scan ------------------------------------------------------


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



# xp-ehh analysis (cross population) --------------------------------------

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
