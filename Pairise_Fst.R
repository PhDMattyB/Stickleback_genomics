##############################
## Fst b/w ecotypes
##
## Matt Brachmann (PhDMattyB)
##
## 14.05.2024
##
##############################


setwd('~/Parsons_Postdoc/Stickleback_Genomic/vcf_filter/')

library(ape)
library(adegenet)
library(pegas)
library(hierfstat)
library(poppr)
library(tidyverse)
library(vcfR)
library(StAMPP)
# 
# data(jaguar)
# Fst(jaguar)

vcf = read.vcfR('stickleback_filtered_vcf.vcf', 
                verbose = FALSE)
genlight = vcfR2genlight(vcf)

meta_df = read_csv('Whole_genome_IndividualID.csv') %>% 
  rename(none = 2, 
         pops = 3)

genlight@pop = as.factor(meta_df$pops)
# stickle_data@pop = as.factor(meta_df$pops)
# stickle_data@pop

stickle_pops = popsub(genlight, 
                  sublist = c('ASHNC', 
                              'ASHNW', 
                              'MYVC', 
                              'MYVW', 
                              'SKRC', 
                              'SKRW', 
                              'GTS', 
                              'CSWY'))

fst_cal = stamppFst(geno = stickle_pops, 
          nboots = 100, 
          percent = 95, 
          nclusters = 1)


