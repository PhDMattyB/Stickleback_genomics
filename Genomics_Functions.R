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

Chr_convert = function(data){
  data = mutate(.data = data,
                chr_num = as.factor(case_when(
                  Chromosome == 'chr_I' ~ '1',
                  Chromosome == 'chr_II' ~ '2',
                  Chromosome == 'chr_III' ~ '3',
                  Chromosome == 'chr_IV' ~ '4',
                  Chromosome == 'chr_IX' ~ '5',
                  Chromosome == 'chr_M' ~ '6',
                  Chromosome == 'chr_Un' ~ '7',
                  Chromosome == 'chr_V' ~ '8',
                  Chromosome == 'chr_VI' ~ '9',
                  Chromosome == 'chr_VII' ~ '10',
                  Chromosome == 'chr_VIII' ~ '11',
                  Chromosome == 'chr_X' ~ '12',
                  Chromosome == 'chr_XI' ~ '13',
                  Chromosome == 'chr_XII' ~ '14',
                  Chromosome == 'chr_XIII' ~ '15',
                  Chromosome == 'chr_XIV' ~ '16',
                  Chromosome == 'chr_XIX' ~ '17',
                  Chromosome == 'chr_XV' ~ '18',
                  Chromosome == 'chr_XVI' ~ '19',
                  Chromosome == 'chr_XVII' ~ '20',
                  Chromosome == 'chr_XVIII' ~ '21',
                  Chromosome == 'chr_XX' ~ '22',
                  Chromosome == 'chr_XXI' ~ '23',
                  Chromosome == 'chr_Y' ~ '24')))
  return(data)  
}

Chr_convert2 = function(data){
  data = mutate(.data = data,
                chr_num = as.factor(case_when(
                  CHR == 'chr_I' ~ '1',
                  CHR == 'chr_II' ~ '2',
                  CHR == 'chr_III' ~ '3',
                  CHR == 'chr_IV' ~ '4',
                  CHR == 'chr_IX' ~ '5',
                  CHR == 'chr_M' ~ '6',
                  CHR == 'chr_Un' ~ '7',
                  CHR == 'chr_V' ~ '8',
                  CHR == 'chr_VI' ~ '9',
                  CHR == 'chr_VII' ~ '10',
                  CHR == 'chr_VIII' ~ '11',
                  CHR == 'chr_X' ~ '12',
                  CHR == 'chr_XI' ~ '13',
                  CHR == 'chr_XII' ~ '14',
                  CHR == 'chr_XIII' ~ '15',
                  CHR == 'chr_XIV' ~ '16',
                  CHR == 'chr_XIX' ~ '17',
                  CHR == 'chr_XV' ~ '18',
                  CHR == 'chr_XVI' ~ '19',
                  CHR == 'chr_XVII' ~ '20',
                  CHR == 'chr_XVIII' ~ '21',
                  CHR == 'chr_XX' ~ '22',
                  CHR == 'chr_XXI' ~ '23',
                  CHR == 'chr_Y' ~ '24')))
  return(data)  
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
  # data$CHR = as.character(data$CHR)
  return(data)
  
}

dist_cal = function(data){
  data %>% 
    group_by(CHR) %>% 
    summarise(chr_len = max(POS)) %>% 
    mutate(total = cumsum(chr_len)-chr_len) %>% 
    dplyr::select(-chr_len) %>% 
    left_join(data, ., by = c('CHR'='CHR')) %>%
    arrange(CHR, 
            POS) %>% 
    mutate(BPcum = POS + total) 
}


axis_df = function(data){
  data %>% 
    group_by(CHR) %>% 
    summarize(center=(max(BPcum) + min(BPcum))/2 )  
}



Fst_manhattan = function(non_outs, 
                         outs, 
                         axisdf,
                         xval, 
                         yval, 
                         chr, 
                         out_col = '', 
                         plot_letter = ''){
  ggplot(non_outs, 
         aes(x = {{xval}}, 
             y = {{yval}}))+
    # plot the non outliers in grey
    geom_point(aes(color = as.factor(chr)), 
               alpha = 0.8, 
               size = 1.3)+
    ## alternate colors per chromosome
    scale_color_manual(values = rep(c("grey", "dimgrey"), 39))+
    ## plot the outliers on top of everything
    ## currently digging this hot pink colour
    geom_point(data = outs,
               col = out_col,
               alpha=0.8, 
               size=1.3)+
    scale_x_continuous(label = axisdf$CHR, 
                       breaks = axisdf$center)+
    scale_y_continuous(expand = c(0, 0), 
                       limits = c(0,1.0))+
    # geom_hline(yintercept = 0.00043, 
    #            linetype = 2, 
    #            col = 'Black')+
    # ylim(0,1.0)+
    # scale_y_reverse(expand = c(0, 0))+
    # remove space between plot area and x axis
    labs(x = 'Cumulative base pair', 
         y = 'Fst', 
         title = plot_letter)+
    theme(legend.position="none",
          # panel.border = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(), 
          axis.text.x = element_text(size = 9, 
                                     angle = 90), 
          axis.title = element_text(size = 14),
          axis.title.x = element_blank(),
          axis.text.y = element_text(size = 12))
}


LFMM_manhattan = function(non_outs, 
                         outs, 
                         axisdf,
                         xval, 
                         yval, 
                         chr, 
                         out_col = '', 
                         plot_letter = ''){
  ggplot(non_outs, 
         aes(x = {{xval}}, 
             y = {{yval}}))+
    # plot the non outliers in grey
    geom_point(aes(color = as.factor(chr)), 
               alpha = 0.8, 
               size = 1.3)+
    ## alternate colors per chromosome
    scale_color_manual(values = rep(c("grey", "dimgrey"), 39))+
    ## plot the outliers on top of everything
    ## currently digging this hot pink colour
    geom_point(data = outs,
               col = out_col,
               alpha=0.8, 
               size=1.3)+
    scale_x_continuous(label = axisdf$CHR, 
                       breaks = axisdf$center)+
    # scale_y_continuous(expand = c(0, 0), 
    #                    limits = c(0,1.0))+
    # geom_hline(yintercept = 0.00043, 
    #            linetype = 2, 
    #            col = 'Black')+
    # ylim(0,1.0)+
    scale_y_reverse(expand = c(0, 0), 
                    limits = c(0,1.0))+
    # remove space between plot area and x axis
    labs(x = 'Cumulative base pair', 
         y = 'qvalue', 
         title = plot_letter)+
    theme(legend.position="none",
          # panel.border = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(), 
          axis.text.x = element_text(size = 9, 
                                     angle = 90), 
          axis.title = element_text(size = 14),
          axis.title.x = element_blank(),
          axis.text.y = element_text(size = 12))
}

# this function takes a .gff or .gff3 file and creates a tibble
# arguments:
# filename -  the name of the gff file
# lead_skip - this argument specifies the
# chromosome - the chromosome to subset from the full gff, default is 'all' and no subset is performed.
read_gff = function(filename, lead_skip = 9, chromosome = "all"){
  #CAM - personally I would pass a vector of names to the read_tsv as opposed to renaming in the pipe (this could be an argument too).
  colnames = c("chr", "source", "feature", "start", "end", "score", "strand", "frame", "attribute" )
  #read in the raw file with readr
  gff = read_tsv(filename, col_names = colnames, skip = lead_skip) %>%
    filter(feature %in% c('gene', 'CDS')) %>%
    arrange(start, end) %>% 
    mutate(mid = start + (end-start)/2) ## arrange each gene by its start and end points on each chromosome
  #subset if user changed the chromosome argument
  if(chromosome != "all"){
    sub_gff = gff %>% filter(chr == chromosome)
    return(sub_gff)
  }
  #if the default chr arg was not changed, return full gff
  return(gff)
}


# Sliding window functions ------------------------------------------------


Mb_Conversion = function(data){
  data %>% 
    group_by(AC_CHR) %>% 
    mutate(win_mid_mb = win_mid/1000000)
}

SW_dist_cal = function(data){
  data %>% 
    group_by(CHR) %>% 
    summarise(chr_len = max(win_mid)) %>% 
    mutate(total = cumsum(chr_len)-chr_len) %>% 
    dplyr::select(-chr_len) %>% 
    left_join(data, ., by = c('CHR'='CHR')) %>%
    arrange(CHR, 
            win_mid) %>% 
    mutate(BPcum = win_mid + total) 
}


SW_top_5_outliers = function(data){
 
   data[data$FST_mean > quantile(data$FST_mean, 
                                prob = 1-5/100),]
}

SW_top_0.5_outliers = function(data){
  
  data[data$FST_mean > quantile(data$FST_mean, 
                                prob = 1-0.5/100),]
}


# Create tped -------------------------------------------------------------

Create_tped = function(ped, map){
  ## Obtaining a vector of names for each individual
  message('Creating vector of individual names')
  indiv_names = ped %>% 
    # filter(FamilyID != 'GDL') %>%
    select(IndividualID) %>% 
    t() %>% 
    as_tibble() %>% 
    row_to_names(row_number = 1) %>% 
    names()
  
  message('Merging ped and map files')
  ped %>% 
    # filter(FamilyID != 'GDL') %>%
    select(contains('AX-')) %>% 
    t() %>% 
    as_tibble() %>% 
    rename_all(funs(c(indiv_names))) %>% 
    bind_cols(map) %>% 
    select(Chromosome, 
           MarkerID, 
           Genetic_dist, 
           Physical_dist, 
           everything())
}
