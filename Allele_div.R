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

write.vcf(x = chrxxi_data, 
          'chrxxi_vcf.vcf')

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

chrxxi_dart = dartR.base::gl.read.vcf('chrxxi_vcf.vcf')

# chr_XXI_SNPs$POS = as.numeric(chr_XXI_SNPs$POS)
# 
# chr_xxi_invert = chr_XXI_SNPs %>% 
#   dplyr::filter(POS >= 9963830,
#                 POS <= 11574024)
# 
# chrxxi_loci = vcfR2loci(chrxxi_data)
# 
# chrxxi_genind = vcfR2genind(chrxxi_data)
# chrxxi_genlight = vcfR2genlight(chrxxi_data)
# 
# glPca(chrxxi_genlight)
# 
# chrxxi_genlight@position
# 
# gl <- dartR.base::gl.compliance.check(chrxxi_genlight)
# gl <- dartR.base::gl.recalc.metrics(chrxxi_genlight)

# chrxxi_genlight@other$loc.metrics = chrxxi_genlight$n.loc
# 
# num_loci_genlight <- ncol(chrxxi_genlight)
# num_loci_metrics <- nrow(chrxxi_genlight$other$loc.metrics)

chrxxi_dart@other$loc.metrics$position = chrxxi_dart@position


filtered_gl_data = dartR.base::gl.filter.locmetric(
  x = chrxxi_dart,
  metric = "position", # The name of the metric you added
  lower = 9963830,
  upper = 11574024,
  keep = "within",
  verbose = 3 # For detailed output
)

pca_inversion = adegenet::glPca(filtered_gl_data, 
                      nf = 2)

pca_inversion$loadings
pca_inversion$scores
pca_inversion$eig

inversion_pca_load = pca_inversion$loadings %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  as.tibble() %>% 
  rename(SNP = rowname)

inversion_pca_score = pca_inversion$scores %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  as.tibble() %>% 
  rename(Individual = rowname)

ggplot(data = inversion_pca_score, 
       aes(x = PC1, 
           y = PC2))+
  geom_point()

