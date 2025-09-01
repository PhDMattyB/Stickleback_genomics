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
library(ChromHeatMap)

vcf = read.vcfR('stickleback_filtered_vcf.vcf')

head(vcf)

chrxxi_data = vcf[getCHROM(vcf) == 'chr_XXI']

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


dna <- ape::read.dna('GCF_016920845.1_GAculeatus_UGA_version5_genomic.fna', 
                     format = "fasta")

gff <- read.table('genomic.gff', sep="\t", quote="")

chromR_xxi = create.chromR(name = 'chrxxi', 
                           vcf = chrxxi_data, 
                           seq = dna, 
                           ann = gff, 
                           verbose = T)

plot(chromR_xxi)

chromoqc(chromR_xxi, 
         dp.alpha = 66)

proc_chromR_xxi = proc.chromR(chromR_xxi, verbose = TRUE)

plot(proc_chromR_xxi)
chromoqc(proc_chromR_xxi)
