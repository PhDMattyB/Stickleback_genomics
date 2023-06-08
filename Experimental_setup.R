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
           Rep, 
           Temperature)%>% 
  rename(Pop = Population, 
         Temp = Temperature) %>% 
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
           Rep, 
           Temperature)%>% 
  rename(Pop = Population, 
         Temp = Temperature) %>% 
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
           Rep, 
           Temp)%>% 
  summarize(avg_weight_final = mean(Weight), 
            avg_length_final = mean(Length))

low_food_warm_final = Final_data %>% 
  # filter(Treatment == 'Low Food', 
  #        Population == 'ACAC') %>% 
  filter(Treatment == 'Low Food', 
         Temp == '18') %>% 
  group_by(Pop, 
           Rep, 
           Temp)%>% 
  summarize(avg_weight_final = mean(Weight), 
            avg_length_final = mean(Length))


# Warm Data Comparison ----------------------------------------------------


high_warm_results = inner_join(high_food_warm, 
           high_food_warm_final, 
           by = c('Pop', 
                  'Rep', 
                  'Temp')) %>% 
  select(-'high_food_amount') %>% 
  mutate(Weight_dif = avg_weight_final - avg_weight, 
         Length_growth = avg_length_final - avg_length) %>% 
  select(Pop, 
         Rep,
         Temp,
         Weight_dif, 
         Length_growth)


low_warm_results = inner_join(low_food_warm, 
           low_food_warm_final, 
           by = c('Pop', 
                  'Rep', 
                  'Temp')) %>% 
  select(-'low_food_amount') %>% 
  mutate(Weight_dif = avg_weight_final - avg_weight, 
         Length_growth = avg_length_final - avg_length) %>% 
  select(Pop, 
         Rep, 
         Temp,
         Weight_dif, 
         Length_growth)


# Cold side ---------------------------------------------------------------


low_food_cold = treatment_data %>% 
  # filter(Treatment == 'Low Food', 
  #        Population == 'ACAC') %>% 
  filter(Treatment == 'Low Food', 
         Temperature == '12') %>% 
  group_by(Population, 
           Rep, 
           Temperature)%>% 
  rename(Pop = Population, 
         Temp = Temperature) %>% 
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
           Rep, 
           Temperature)%>% 
  rename(Pop = Population, 
         Temp = Temperature) %>% 
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
           Rep, 
           Temp)%>% 
  summarize(avg_weight_final = mean(Weight), 
            avg_length_final = mean(Length))

low_food_cold_final = Final_data %>%
  # filter(Treatment == 'Low Food',
  #        Population == 'ACAC') %>%
  filter(Treatment == 'Low Food',
         Temp == '12') %>%
  group_by(Pop,
           Rep, 
           Temp)%>%
  summarize(avg_weight_final = mean(Weight),
            avg_length_final = mean(Length))



# Cold data compare --------------------------------------------------------

high_cold_results = inner_join(high_food_cold, 
                               high_food_cold_final, 
                               by = c('Pop', 
                                      'Rep', 
                                      'Temp')) %>% 
  select(-'high_food_amount') %>% 
  mutate(Weight_dif = avg_weight_final - avg_weight, 
         Length_growth = avg_length_final - avg_length) %>% 
  select(Pop, 
         Rep, 
         Temp,
         Weight_dif, 
         Length_growth)


low_cold_results = inner_join(low_food_cold,
                              low_food_cold_final,
                              by = c('Pop',
                                     'Rep', 
                                     'Temp')) %>%
  select(-'low_food_amount') %>%
  mutate(Weight_dif = avg_weight_final - avg_weight,
         Length_growth = avg_length_final - avg_length) %>%
  select(Pop,
         Rep,
         Temp,
         Weight_dif,
         Length_growth)


# Cold warm comparison ----------------------------------------------------

high_warm_results = high_warm_results %>% 
  rename(weight_dif_warm = Weight_dif, 
         Length_dif_warm = Length_growth)

high_cold_results = high_cold_results %>% 
  rename(weight_dif_cold = Weight_dif, 
         Length_dif_cold = Length_growth)

high_compare_results = inner_join(high_warm_results, 
           high_cold_results, 
           by = c('Pop', 
                  'Rep')) %>% 
  mutate(high_warm_weight_diff = weight_dif_warm - weight_dif_cold, 
         high_warm_length_diff = Length_dif_warm - Length_dif_cold) %>% 
  dplyr::select(Pop, 
                Rep,
                high_warm_weight_diff, 
                high_warm_length_diff)





low_warm_results = low_warm_results %>% 
  rename(weight_dif_warm = Weight_dif, 
         Length_dif_warm = Length_growth)

low_cold_results = low_cold_results %>% 
  rename(weight_dif_cold = Weight_dif, 
         Length_dif_cold = Length_growth)

low_compare_results = inner_join(low_warm_results, 
           low_cold_results, 
           by = c('Pop', 
                  'Rep')) %>% 
  mutate(low_warm_weight_diff = weight_dif_warm - weight_dif_cold, 
         low_warm_length_diff = Length_dif_warm - Length_dif_cold) %>% 
  dplyr::select(Pop, 
                Rep,
                low_warm_weight_diff, 
                low_warm_length_diff)


high_compare_results

s