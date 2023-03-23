##############################
## Experimental set up
##
## Matt Brachmann (PhDMattyB)
##
## 23.03.2023 
##
##############################


setwd('~/Parsons_Postdoc/Experiment1/')

library(tidyverse)

treatment_data = read_csv('Experiment1_setup_data.csv')


treatment_data %>% 
  group_by(Treatment, 
           Rep) %>% 
  summarize(avg_weight = mean(Weight), 
            avg_length = mean(Length)) %>% 
  mutate(low_food_amount = avg_weight*0.02, 
         high_food_amount = avg_weight*0.3)
