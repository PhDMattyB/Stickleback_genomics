##############################
## Fish tagging fiasco
##
## Matt Brachmann (PhDMattyB)
##
## 31.07.2023
##
##############################

library(tidyverse)

setwd('~/Parsons_Postdoc/')

fiasco_data = read_csv('Scarred_fish_data.csv')

mean_vals = fiasco_data %>% 
  group_by(Cross) %>% 
  mutate(mean_weight = mean(weight_gm), 
         mean_length = mean(length_mm)) %>%
  summarize(mean_weight = mean(weight_gm), 
            mean_length = mean(length_mm)) 
