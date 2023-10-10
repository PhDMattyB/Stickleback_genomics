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



# Haplotype structure around selection target -----------------------------
# find the highest hit
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
# find the highest hit
hit = warm_cold %>% 
  arrange(desc(LOGPVALUE)) %>% 
  top_n(1)

# get SNP position
x = hit$POSITION

marker_id_warm = which(warm@positions == x)
marker_id_cold = which(cold@positions == x)

warm_furcation = calc_furcation(warm, 
                                mrk = marker_id_warm)

cold_furcation = calc_furcation(cold, 
                                mrk = marker_id_cold )


plot(warm_furcation, 
     xlim = c(58695, 17420697))
plot(cold_furcation, 
     xlim = c(58694, 17420697))

warm_haplen = calc_haplen(warm_furcation)
cold_haplen = calc_haplen(cold_furcation)

plot(warm_haplen)
plot(cold_haplen)


# Warm vs Cold sweep detect (chr 21) --------------------------------------

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

## plot the ihs statistic
ggplot(warm2_ihs$ihs, 
       aes(POSITION, 
           IHS))+
  geom_point()

ggplot(cold2_ihs$ihs, 
       aes(POSITION, 
           IHS))+
  geom_point()

## plot the pvalues
ggplot(warm2_ihs$ihs, 
       aes(POSITION, 
           LOGPVALUE))+
  geom_point()

ggplot(cold2_ihs$ihs, 
       aes(POSITION, 
           LOGPVALUE))+
  geom_point()



# xp-ehh analysis (cross population) 

warm2_cold2 = ies2xpehh(warm2_scan, 
                      cold2_scan, 
                      popname1 = 'warm2', 
                      popname2 = 'cold2', 
                      include_freq = T)

# plot

ggplot(warm2_cold2, 
       aes(POSITION, 
           XPEHH_warm2_cold2))+
  geom_point()

ggplot(warm2_cold2, 
       aes(POSITION, 
           LOGPVALUE)) + geom_point()



# Haplotype structure around selection target 
# find the highest hit
hit = warm2_cold2 %>% 
  arrange(desc(LOGPVALUE)) %>% 
  top_n(1)

# get SNP position
x = hit$POSITION

marker_id_warm2 = which(warm2@positions == x)
marker_id_cold2 = which(cold2@positions == x)

warm2_furcation = calc_furcation(warm2, 
                                mrk = marker_id_warm2)

cold2_furcation = calc_furcation(cold2, 
                                mrk = marker_id_cold2 )


plot(warm2_furcation, 
     xlim = c(58695, 17420697))
plot(cold2_furcation, 
     xlim = c(58694, 17420697))

warm2_haplen = calc_haplen(warm2_furcation)
cold2_haplen = calc_haplen(cold2_furcation)

plot(warm2_haplen)
plot(cold2_haplen)


