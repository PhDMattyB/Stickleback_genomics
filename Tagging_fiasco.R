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

## The fucking 'problem' cross is bigger than the 'normals'

weight_anova = aov(weight_gm ~ Cross, 
                   data = fiasco_data)

summary(weight_anova)

length_anova = aov(length_mm ~ Cross, 
                   data = fiasco_data)

summary(length_anova)

