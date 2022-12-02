##############################
## Genomic functions misc. 
##
## Matt Brachmann (PhDMattyB)
##
## 2022-12-02
##
##############################

library(tidyverse)

theme_set(theme_bw())

Fst_manhatan_format = function(Fst_data, Fst_outliers){
  Neutral_snps = anti_join(Fst_data, 
                           Fst_outliers)
  ## Need to remove outliers from this dataset
  labs1 = rep('Neutral', 
              nrow(Neutral_snps)) %>% 
    as_tibble()
  
  Fst_neutral = bind_cols(Neutral_snps, 
                          labs1)
  
  
  labs2 = rep('Outlier', 
              nrow(Fst_outliers)) %>% 
    as_tibble()
  
  Fst_outliers = bind_cols(Fst_outliers, 
                           labs2)
  
  
  Fst_clean = bind_rows(Fst_neutral, 
                        Fst_outliers)
  
  
}


stickle_CHR_reorder = function(data){
  data$CHR <- factor(data$CHR, 
                                   levels=c("chr_I",
                                            "chr_II", 
                                            "chr_III", 
                                            'chr_IV', 
                                            'chr_V', 
                                            'chr_VI', 
                                            'chr_VII', 
                                            'chr_VIII', 
                                            'chr_IX', 
                                            'chr_X', 
                                            'chr_XI', 
                                            'chr_XII', 
                                            'chr_XIII', 
                                            'chr_XIV', 
                                            'chr_XV', 
                                            'chr_XVI', 
                                            'chr_XVII', 
                                            'chr_XVIII', 
                                            'chr_XIX', 
                                            'chr_XX', 
                                            'chr_XXI', 
                                            'chr_Y', 
                                            'chr_M', 
                                            'chr_Un'))
  
}
