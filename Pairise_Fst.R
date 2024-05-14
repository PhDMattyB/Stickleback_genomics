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
library(tidyverse)

data(jaguar)
Fst(jaguar)


# read.vcf('stickleback_filtered_vcf.vcf', 
#          which.loci = 1:173000)

# read.genepop('stickleback_filtered_genepop.gen')

stickle_data = read.structure('stickleback_filtered_structure.stru')
meta_df = read_csv('Whole_genome_IndividualID.csv') %>% 
  rename(none = 2, 
         pops = 3)

stickle_data@pop = as.factor(meta_df$pops)
# stickle_data@pop

pegas_df = alleles2loci(stickle_data)

Fst(pegas_df, 
    pop = meta_df$pops)
