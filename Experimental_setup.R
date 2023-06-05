##############################
## Experimental set up
##
## Matt Brachmann (PhDMattyB)
##
## 23.03.2023 
##
##############################


setwd('~/Parsons_Postdoc/Experiment1/')

# setwd('~/Parsons_Postdoc/Experiment1/')

library(tidyverse)

treatment_data = read_csv('Experiment1_setup_data.csv')

##
# Warm side Initial---------------------------------------------------------------


low_food_warm = treatment_data %>% 
  # filter(Treatment == 'Low Food', 
  #        Population == 'ACAC') %>% 
  filter(Treatment == 'Low Food', 
         Temperature == '18') %>% 
  group_by(Population, 
           Rep)%>% 
  rename(Pop = Population) %>% 
  summarize(avg_weight = mean(Weight), 
            avg_length = mean(Length)) %>% 
  mutate(low_food_amount = avg_weight*0.02)

treatment_data %>% 
  filter(Treatment == 'Low Food', 
         Temperature == '18') %>% 
  group_by(Population, 
           Rep) %>% 
  summarize(Number = n())


high_food_warm = treatment_data %>% 
  # filter(Treatment == 'High Food', 
  #        Population == 'ACAC') %>% 
  filter(Treatment == 'High Food', 
         Temperature == '18') %>% 
  group_by(Population, 
           Rep)%>% 
  rename(Pop = Population) %>% 
  summarize(avg_weight = mean(Weight), 
            avg_length = mean(Length)) %>% 
  mutate(high_food_amount = avg_weight*0.3)

treatment_data %>% 
  filter(Treatment == 'High Food', 
         Temperature == '18') %>% 
  group_by(Population, 
           Rep) %>% 
  summarize(Number = n())


# Warm side final ---------------------------------------------------------


Final_data = read_csv('Experiment1_Sampling_Good.csv')

high_food_warm_final = Final_data %>% 
  # filter(Treatment == 'High Food', 
  #        Population == 'ACAC') %>% 
  filter(Treatment == 'High Food', 
         Temp == '18') %>% 
  group_by(Pop, 
           Rep)%>% 
  summarize(avg_weight_final = mean(Weight), 
            avg_length_final = mean(Length))

low_food_warm_final = Final_data %>% 
  # filter(Treatment == 'Low Food', 
  #        Population == 'ACAC') %>% 
  filter(Treatment == 'Low Food', 
         Temp == '18') %>% 
  group_by(Pop, 
           Rep)%>% 
  summarize(avg_weight_final = mean(Weight), 
            avg_length_final = mean(Length))


# Warm Data Comparison ----------------------------------------------------


high_warm_results = inner_join(high_food_warm, 
           high_food_warm_final, 
           by = c('Pop', 
                  'Rep')) %>% 
  select(-'high_food_amount') %>% 
  mutate(Weight_dif = avg_weight_final - avg_weight, 
         Length_growth = avg_length_final - avg_length) %>% 
  select(Pop, 
         Rep, 
         Weight_dif, 
         Length_growth)


low_warm_results = inner_join(low_food_warm, 
           low_food_warm_final, 
           by = c('Pop', 
                  'Rep')) %>% 
  select(-'low_food_amount') %>% 
  mutate(Weight_dif = avg_weight_final - avg_weight, 
         Length_growth = avg_length_final - avg_length) %>% 
  select(Pop, 
         Rep, 
         Weight_dif, 
         Length_growth)


# Cold side ---------------------------------------------------------------


low_food_cold = treatment_data %>% 
  # filter(Treatment == 'Low Food', 
  #        Population == 'ACAC') %>% 
  filter(Treatment == 'Low Food', 
         Temperature == '12') %>% 
  group_by(Population, 
           Rep)%>% 
  rename(Pop = Population) %>% 
  summarize(avg_weight = mean(Weight), 
            avg_length = mean(Length)) %>% 
  mutate(low_food_amount = avg_weight*0.02)

treatment_data %>% 
  filter(Treatment == 'Low Food', 
         Temperature == '12') %>% 
  group_by(Population, 
           Rep) %>% 
  summarize(Number = n())


high_food_cold = treatment_data %>% 
  # filter(Treatment == 'High Food', 
  #        Population == 'ACAC') %>% 
  filter(Treatment == 'High Food', 
         Temperature == '12') %>% 
  group_by(Population, 
           Rep)%>% 
  rename(Pop = Population) %>% 
  summarize(avg_weight = mean(Weight), 
            avg_length = mean(Length)) %>% 
  mutate(high_food_amount = avg_weight*0.3)

treatment_data %>% 
  filter(Treatment == 'High Food', 
         Temperature == '12') %>% 
  group_by(Population, 
           Rep) %>% 
  summarize(Number = n())



# Cold side final ---------------------------------------------------------

Final_data = read_csv('Experiment1_Sampling_Good.csv')

high_food_cold_final = Final_data %>% 
  # filter(Treatment == 'High Food', 
  #        Population == 'ACAC') %>% 
  filter(Treatment == 'High Food', 
         Temp == '12') %>% 
  group_by(Pop, 
           Rep)%>% 
  summarize(avg_weight_final = mean(Weight), 
            avg_length_final = mean(Length))

low_food_cold_final = Final_data %>%
  # filter(Treatment == 'Low Food',
  #        Population == 'ACAC') %>%
  filter(Treatment == 'Low Food',
         Temp == '18') %>%
  group_by(Pop,
           Rep)%>%
  summarize(avg_weight_final = mean(Weight),
            avg_length_final = mean(Length))



# Cold data compare --------------------------------------------------------

high_cold_results = inner_join(high_food_cold, 
                               high_food_cold_final, 
                               by = c('Pop', 
                                      'Rep')) %>% 
  select(-'high_food_amount') %>% 
  mutate(Weight_dif = avg_weight_final - avg_weight, 
         Length_growth = avg_length_final - avg_length) %>% 
  select(Pop, 
         Rep, 
         Weight_dif, 
         Length_growth)


low_cold_results = inner_join(low_food_cold,
                              low_food_cold_final,
                              by = c('Pop',
                                     'Rep')) %>%
  select(-'low_food_amount') %>%
  mutate(Weight_dif = avg_weight_final - avg_weight,
         Length_growth = avg_length_final - avg_length) %>%
  select(Pop,
         Rep,
         Weight_dif,
         Length_growth)


# Cold warm comparison ----------------------------------------------------

high_warm_results

high_cold_results
