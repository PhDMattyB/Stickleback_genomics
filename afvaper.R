library(dplyr)
library(data.table)
library(vcfR)
library(afvaper)
library(gtools)
library(readr)


# setwd('~/Documents/afvaper')
# setwd('~/Documents/afvaper/no_chr23/')
setwd('~/Parsons_Postdoc/Stickleback_Genomic/afvaper/')
##
# afvaper single chromosome -----------------------------------------------
## single chromosome
chr_vcf = read.vcfR('stickle_filtered_24.vcf')


## individual ids got duplicated in the per chr vcfs somehow
# chr1_vcf@gt

popmap = read_tsv('Stickleback_afvaper_popmap_dup.txt', 
                  col_names = F) %>% 
  as.data.frame()

input_vectors = list(pair1 = c("Warm1", "Cold1"),
                     pair2 = c("Warm2", "Cold2"),
                     pair3 = c("Warm3", "Cold3"),
                     pair4 = c("Warm4", "Cold4"))


# Set our window size
window_snps = 50

# Calculate Allele Frequency Change Vector Matrices
chr_af = calc_AF_vectors(vcf = chr_vcf,
                           window_size = window_snps,
                           popmap = popmap,
                           vectors = input_vectors,
                           n_cores = 1,
                           data_type = "vcf")




# How many permutations to run

null_perm_N = 331

# Calculate Allele Frequency Change Vector Matrices
null_input = calc_AF_vectors(vcf = chr_vcf,
                             window_size = window_snps,
                             popmap = popmap,
                             vectors = input_vectors,
                             n_cores = 1,
                             null_perms = null_perm_N,
                             data_type = "vcf")


# afvaper multiple chromosomes --------------------------------------------

chr_perms = read_tsv('Stickleback_Chromosome_sizes_fixed.txt') %>% 
  slice(-6)
# chr_perms = read_tsv('Stickleback_Chromosome_sizes_100000.txt') %>% 
#   slice(-6)
window_snps = 10
# window_snps = 50
# window_snps = 100
# window_snps = 200
popmap = read_tsv('Stickleback_afvaper_popmap_dup.txt', 
                  col_names = F) %>% 
  as.data.frame()

input_vectors = list(pair1 = c("Warm1", "Cold1"),
                     pair2 = c("Warm2", "Cold2"),
                     pair3 = c("Warm3", "Cold3"),
                     pair4 = c("Warm4", "Cold4"))

## check the number of chrs with the lapply numbers
all_chr_res <- lapply(1:23,function(i){
  
  # First read in the VCF
  list_vcf_files = list.files(path = "~/Parsons_Postdoc/Stickleback_Genomic/afvaper/afvaper_chrs/")
  # list_vcf_files = list.files(path = "~/Documents/afvaper/no_chr23/")
  list_vcf_files = mixedsort(list_vcf_files)
  chr_vcf <- read.vcfR(list_vcf_files[i])
  
  # Calculate AFV
  chr_AF_input <- calc_AF_vectors(vcf = chr_vcf,
                                  window_size = window_snps,
                                  popmap = popmap,
                                  vectors = input_vectors,
                                  n_cores = 4)
  
  # Calculate null AFV
  chr_null_input <- calc_AF_vectors(vcf = chr_vcf,
                                    window_size = window_snps,
                                    popmap = popmap,
                                    vectors = input_vectors,
                                    n_cores = 4,
                                    null_perms = chr_perms$perms[i])
  
  ## We could save these to some temporary file, e.g.
  # saveRDS(list(chr_AF_input,chr_null_input),paste0("chr",i,"_AFV.rds"))
  
  # Return our results
  return(list(chr_AF_input,chr_null_input))
})

 # To fetch all of  chr AFV, we take the first element of each list element
# Note: the merge_eigen_res() func is the same as unlist(,recursive=F)
AF_input <- merge_eigen_res(lapply(all_chr_res,'[[',1))

# All null, we take the second element of each list element
null_input <- merge_eigen_res(lapply(all_chr_res,'[[',2))

# We now have our whole (3 chr) genome's worth of AFV matrices and null matrices in a single input
c(head(names(AF_input)),tail(names(AF_input)))

## eigenvector analysis!!
# Perform eigen analysis over allele freq matrix
eigen_res <- lapply(AF_input,eigen_analyse_vectors)

head(names(eigen_res))

eigen_res[[1]]$eigenvals
eigen_res[[1]]$eigenvecs

## null cutoffs
# Get cutoffs for 95%, 99% and 99.9%
null_cutoffs <- find_null_cutoff(null_input,cutoffs = c(0.95,0.99,0.999))
null_cutoffs


# Calculate p-vals
pvals <- eigen_pvals(eigen_res,null_input)
head(pvals)



# Plot the raw eigenvalues, and visualise the cutoff of 99%
all_plots <- eigenval_plot(eigen_res,cutoffs = null_cutoffs[,"99%"])

# Show the plots for eigenvalue 1
eig1_vec__fig = all_plots[[1]]

eig1_vec__fig+
  theme_bw()+
  theme(title = element_blank())



# Plot empirical p-values, -log10(p) of 2 ~ p=0.01, 3 ~ p=0.001 etc.
all_plots_p <- eigenval_plot(eigen_res,
                             null_vectors = null_input,
                             plot.pvalues = T)
all_plots_p[[1]]

### PLOT SPECIFIC CHROMOSOMES!!
# Plot empirical p-values, -log10(p) of 2 ~ p=0.01, 3 ~ p=0.001 etc.
chr_windows <- grep("23",names(eigen_res))
all_plots_p_chr <- eigenval_plot(eigen_res[chr_windows],
                                  null_vectors = null_input,
                                  plot.pvalues = T)

all_plots_p_chr1[[1]]



# pull out significant data per chr ---------------------------------------

# Recall the use of find_null_cutoffs() to fetch a matrix of cutoffs...
# null_cutoffs

# Find significant windows above 99.9% null permutation
significant_windows <- signif_eigen_windows(eigen_res,
                                            null_cutoffs[,"99%"])

# Display 'outliers'
significant_windows



#Summarise parallel evolution in windows that are significant on eigenvector 1
eig1_parallel <- summarise_window_parallelism(window_id = significant_windows[[1]],
                                              eigen_res = eigen_res,
                                              loading_cutoff = 0.3,
                                              eigenvector = 1)
# Show results
head(eig1_parallel)

eig2_parallel <- summarise_window_parallelism(window_id = significant_windows[[2]],
                                              eigen_res = eigen_res,
                                              loading_cutoff = 0.3,
                                              eigenvector = 2)

# Show results
head(eig2_parallel)


#Fetch an A matrix
A_mat <- eigen_res[[1]]$A_matrix
head(A_mat)

eigen_res[[23]]$A_matrix

to_plot <- data.frame(snp=rownames(A_mat),
                      eig1_score=A_mat[,1])

to_plot <- to_plot %>% tidyr::separate("snp",into=c("chr","pos"),sep="_")
to_plot$pos <- as.integer(to_plot$pos)

ggplot(to_plot,aes(x=pos,y=abs(eig1_score)))+
  geom_point()+
  labs(y="Eig1 Score",x="Pos (bp)")
