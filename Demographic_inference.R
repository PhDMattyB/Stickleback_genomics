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



# Per population folded SFS -----------------------------------------------


setwd('~/Parsons_Postdoc/Stickleback_Genomic/vcf_filter/')

MYV = read.vcf('MYV_only.vcf', 
                        from = 1, 
                        to = 173444)

MYV_sfs = site.spectrum(MYV, 
                         folded = T)

plot(MYV_sfs)

MYV_sfs

ASHN = read.vcf('ASHN_only.vcf', 
               from = 1, 
               to = 173444)

ASHN_sfs = site.spectrum(ASHN, 
                        folded = T)


SKR = read.vcf('SKR_only.vcf', 
               from = 1, 
               to = 173444)

SKR_sfs = site.spectrum(SKR, 
                        folded = T)

CSWY = read.vcf('CSWY_only.vcf', 
               from = 1, 
               to = 173444)

CSWY_sfs = site.spectrum(CSWY, 
                        folded = T)

GTS = read.vcf('GTS_only.vcf', 
               from = 1, 
               to = 173444)

GTS_sfs = site.spectrum(GTS, 
                        folded = T)



# Stairway plot 2 results -------------------------------------------------

setwd('~/Parsons_Postdoc/Stickleback_Genomic/stairwayplot2/')


MYV_data = read_tsv('MYV stairway.final.summary') %>% 
  mutate(year_conv = year/1000)
SKR_data = read_tsv('SKR stairway.final.summary')
ASHN_data = read_tsv('ASHN stairway.final.summary')
CSWY_data = read_tsv('CSWY stairway.final.summary')
GTS_data = read_tsv('GTS stairway.final.summary')


View(MYV_data)


ggplot(data = MYV_data, 
       aes(x = year_conv,
           y = Ne_median))+
  geom_line()


ggplot()+
  geom_line(data = MYV_data, 
            aes(x = year, 
                y = Ne_median))+
  geom_line(data = SKR_data, 
            aes(x = year, 
                y = Ne_median), 
            col = 'red')+
  scale_x_reverse()



# nucleotide diversity ----------------------------------------------------

warm_stickle = read.vcf('stickleback_warm_pops.vcf', 
                        from = 1, 
                        to = 173444)

nuc.div(warm_stickle)

