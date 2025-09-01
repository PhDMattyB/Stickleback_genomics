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
# library(ChromHeatMap)

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



# dna <- ape::read.dna('GCF_016920845.1_GAculeatus_UGA_version5_genomic.fna', 
#                      format = "fasta")
# 
# gff <- read.table('genomic.gff', sep="\t", quote="")
# 
# chromR_xxi = create.chromR(name = 'chrxxi', 
#                            vcf = chrxxi_data, 
#                            seq = dna, 
#                            ann = gff, 
#                            verbose = T)
# 
# plot(chromR_xxi)
# 
# chromoqc(chromR_xxi, 
#          dp.alpha = 66)
# 
# proc_chromR_xxi = proc.chromR(chromR_xxi, 
#                               win.size = 50,
#                               verbose = TRUE)
# 
# plot(proc_chromR_xxi)
# chromoqc(proc_chromR_xxi)
# 
# head(proc_chromR_xxi@var.info)
# head(proc_chromR_xxi@win.info)
# 
# 
# chromoqc(proc_chromR_xxi, 
#          xlim=c(9e+06, 12e+06))
# 


# pca on the chr xxi ------------------------------------------------------
library(vegan)
library(dartRverse)

chr_XXI_SNPs$POS = as.numeric(chr_XXI_SNPs$POS)

chr_xxi_invert = chr_XXI_SNPs %>% 
  dplyr::filter(POS >= 9963830,
                POS <= 11574024)

chrxxi_loci = vcfR2loci(chrxxi_data)

chrxxi_genind = vcfR2genind(chrxxi_data)
chrxxi_genlight = vcfR2genlight(chrxxi_data)

glPca(chrxxi_genlight)

chrxxi_genlight@position

gl <- dartR.base::gl.compliance.check(chrxxi_genlight)
gl <- dartR.base::gl.recalc.metrics(chrxxi_genlight)

chrxxi_genlight@other$loc.metrics = chrxxi_genlight$n.loc

num_loci_genlight <- ncol(chrxxi_genlight)
num_loci_metrics <- nrow(chrxxi_genlight$other$loc.metrics)

filtered_gl_data <- dartR.base::gl.filter.locmetric(
  x = chrxxi_genlight,
  metric = "position", # The name of the metric you added
  lower = 4101,
  upper = 4600,
  keep = "within",
  verbose = 3 # For detailed output
)

# chrxxi_pcadapt = read.pcadapt('chrxxi_vcf.vcf', type = "vcf")
# 
# 
# chr_xxi_invert_snps = chr_xxi_invert %>% 
#   dplyr::select(-CHROM, 
#                 -POS,
#                 -ID, 
#                 -REF, 
#                 -ALT)
# 
# SNPs_scaled_invert = scale(chr_xxi_invert_snps)
# 
# 
# map_data_invert = read_tsv('chr21_inversion_region_fixed.map', 
#                            col_names = c('CHR',
#                                          'SNP',
#                                          'GPOS',
#                                          'POS'))
# 
# 
# ped_data_invert = read_table('chr21_inversion_region_fixed.ped', 
#                              col_names = c('PopulationID',
#                                            'IndividualID',
#                                            'MaternalID',
#                                            'PaternalID',
#                                            'Sex',
#                                            'Phenotype',
#                                            map_data_invert$SNP))
# 
# 
# test = ped_data_invert %>% 
#   dplyr::select(-PopulationID, 
#                 -IndividualID, 
#                 -MaternalID, 
#                 -PaternalID, 
#                 -Sex, 
#                 -Phenotype)


