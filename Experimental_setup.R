##############################
## Experimental set up
##
## Matt Brachmann (PhDMattyB)
##
## 23.03.2023 
##
##############################


setwd('~/Parsons_Postdoc/Experiment1/')

setwd('~/Parsons_Postdoc/Experiment1/')

library(tidyverse)

treatment_data = read_csv('Experiment1_setup_data.csv')


# Warm side ---------------------------------------------------------------


low_food_amount = treatment_data %>% 
  # filter(Treatment == 'Low Food', 
  #        Population == 'ACAC') %>% 
  filter(Treatment == 'Low Food', 
         Temperature == '18') %>% 
  group_by(Population, 
           Rep)%>% 
  summarize(avg_weight = mean(Weight), 
            avg_length = mean(Length)) %>% 
  mutate(low_food_amount = avg_weight*0.02)

treatment_data %>% 
  filter(Treatment == 'Low Food', 
         Temperature == '18') %>% 
  group_by(Population, 
           Rep) %>% 
  summarize(Number = n())


high_food_amount = treatment_data %>% 
  # filter(Treatment == 'High Food', 
  #        Population == 'ACAC') %>% 
  filter(Treatment == 'High Food', 
         Temperature == '18') %>% 
  group_by(Population, 
           Rep)%>% 
  summarize(avg_weight = mean(Weight), 
            avg_length = mean(Length)) %>% 
  mutate(high_food_amount = avg_weight*0.3)

treatment_data %>% 
  filter(Treatment == 'High Food', 
         Temperature == '18') %>% 
  group_by(Population, 
           Rep) %>% 
  summarize(Number = n())



# Cold side ---------------------------------------------------------------


low_food_amount = treatment_data %>% 
  # filter(Treatment == 'Low Food', 
  #        Population == 'ACAC') %>% 
  filter(Treatment == 'Low Food', 
         Temperature == '12') %>% 
  group_by(Population, 
           Rep)%>% 
  summarize(avg_weight = mean(Weight), 
            avg_length = mean(Length)) %>% 
  mutate(low_food_amount = avg_weight*0.02)

treatment_data %>% 
  filter(Treatment == 'Low Food', 
         Temperature == '12') %>% 
  group_by(Population, 
           Rep) %>% 
  summarize(Number = n())


high_food_amount = treatment_data %>% 
  # filter(Treatment == 'High Food', 
  #        Population == 'ACAC') %>% 
  filter(Treatment == 'High Food', 
         Temperature == '12') %>% 
  group_by(Population, 
           Rep)%>% 
  summarize(avg_weight = mean(Weight), 
            avg_length = mean(Length)) %>% 
  mutate(high_food_amount = avg_weight*0.3)

treatment_data %>% 
  filter(Treatment == 'High Food', 
         Temperature == '12') %>% 
  group_by(Population, 
           Rep) %>% 
  summarize(Number = n())

