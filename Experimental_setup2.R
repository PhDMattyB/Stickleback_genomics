##############################
## Experimental set up
##
## Matt Brachmann (PhDMattyB)
##
## 23.03.2023 
##
##############################


# setwd('')

setwd('~/Parsons_Postdoc/Experiment1/')

library(tidyverse)

treatment_data = read_csv('Experiment1_setup_data.csv') 

# treatment_data %>% 
#   filter(Population == 'ACAC', 
#          Crosses == 126)

# treatment_data %>% 
#   group_by(Treatment,
#           Temperature,
#           Rep) %>% 
#   summarize(n = n(), 
#             avg_length = mean(Length)) %>% 
#   View()
##
# Warm side Initial---------------------------------------------------------------

## quantify condition factor based on that new paper we just read
low_food_warm = treatment_data %>% 
  # filter(Treatment == 'Low Food', 
  #        Population == 'ACAC') %>% 
  filter(Treatment == 'Low Food', 
         Temperature == '18') %>% 
  rename(Pop = Population, 
         Temp = Temperature) %>%
  group_by(Rep,
           Pop,
           Temp) %>%  
  summarize(avg_weight = mean(Weight), 
            avg_length = mean(Length), 
            density = n()) %>% 
  mutate(low_food_amount = avg_weight*0.02,
         low_food_amount_density = avg_weight*0.02/density,
         low_cond_fac = avg_weight/avg_length^(1/3)*100)
# 
# View(low_food_warm)
# 
# treatment_data %>% 
#   filter(Treatment == 'Low Food', 
#          Temperature == '18') %>% 
#   group_by(Population, 
#            Rep) %>% 
#   summarize(Number = n())


high_food_warm = treatment_data %>% 
  # filter(Treatment == 'High Food', 
  #        Population == 'ACAC') %>% 
  filter(Treatment == 'High Food', 
         Temperature == '18') %>%
  rename(Pop = Population, 
         Temp = Temperature) %>%
  group_by(Rep,
           Pop,
           Temp)%>% 
  summarize(avg_weight = mean(Weight), 
            avg_length = mean(Length), 
            density = n()) %>% 
  mutate(high_food_amount = avg_weight*0.3,
         high_food_amount_density = (avg_weight*0.3)/density,
         high_cond_fac = avg_weight/avg_length^(1/3)*100)

# treatment_data %>% 
#   filter(Treatment == 'High Food', 
#          Temperature == '18') %>% 
#   group_by(Population, 
#            Rep) %>% 
#   summarize(Number = n())


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
            avg_length_final = mean(Length), 
            high_cond_fac_final = avg_weight_final/avg_length_final^(1/3)*100)

low_food_warm_final = Final_data %>% 
  # filter(Treatment == 'Low Food', 
  #        Population == 'ACAC') %>% 
  filter(Treatment == 'Low Food', 
         Temp == '18') %>% 
  group_by(Pop, 
           Rep, 
           Temp)%>% 
  summarize(avg_weight_final = mean(Weight), 
            avg_length_final = mean(Length), 
            low_cond_fac_final = avg_weight_final/avg_length_final^(1/3)*100)


# Warm Data Comparison ----------------------------------------------------


high_warm_results = inner_join(high_food_warm, 
           high_food_warm_final, 
           by = c('Pop', 
                  'Rep', 
                  'Temp')) %>% 
  select(-'high_food_amount') %>% 
  mutate(Weight_dif = avg_weight_final - avg_weight, 
         Length_growth = avg_length_final - avg_length, 
         cond_fac_diff = high_cond_fac_final - high_cond_fac) %>% 
  select(Pop, 
         Rep,
         Temp,
         Weight_dif, 
         Length_growth, 
         cond_fac_diff)


low_warm_results = inner_join(low_food_warm, 
           low_food_warm_final, 
           by = c('Pop', 
                  'Rep', 
                  'Temp')) %>% 
  select(-'low_food_amount') %>% 
  mutate(Weight_dif = avg_weight_final - avg_weight, 
         Length_growth = avg_length_final - avg_length, 
         cond_fac_diff = low_cond_fac_final - low_cond_fac) %>% 
  select(Pop, 
         Rep, 
         Temp,
         Weight_dif, 
         Length_growth, 
         cond_fac_diff)


# warm side expectation ---------------------------------------------------

inner_join(low_food_warm, 
           low_food_warm_final, 
           by = c('Pop', 
                  'Rep', 
                  'Temp')) %>% 
  select(-'low_food_amount') %>% 
  mutate(Weight_dif = avg_weight_final - avg_weight, 
         Length_growth = avg_length_final - avg_length, 
         cond_fac_diff = low_cond_fac_final - low_cond_fac)%>% 
  select(Pop, 
         Rep, 
         Temp,
         Weight_dif, 
         Length_growth, 
         cond_fac_diff) %>% 
  ungroup() %>%
  summarize(mean_weight_dif = mean(Weight_dif), 
            mean_length_dif = mean(Length_growth), 
            mean_cond_dif = mean(cond_fac_diff)) %>% 
  mutate(percent_weight_loss = mean_weight_dif*100, 
         percent_length_change = mean_length_dif*100)
  

cold_morph_weight_loss_warm = inner_join(low_food_warm, 
           low_food_warm_final, 
           by = c('Pop', 
                  'Rep', 
                  'Temp')) %>% 
  select(-'low_food_amount') %>% 
  mutate(Weight_dif = avg_weight_final - avg_weight, 
         Length_growth = avg_length_final - avg_length, 
         cond_fac_diff = low_cond_fac_final - low_cond_fac)%>% 
  select(Pop, 
         Rep, 
         Temp,
         Weight_dif, 
         Length_growth, 
         cond_fac_diff) %>% 
  ungroup() %>% 
  filter(Pop %in% c('ACAC', 
                    'MCMC')) %>% 
  summarize(mean_weight_dif = mean(Weight_dif), 
            mean_length_dif = mean(Length_growth), 
            mean_cond_dif = mean(cond_fac_diff)) %>%
  mutate(percent_weight_loss = mean_weight_dif*100, 
         percent_length_change = mean_length_dif*100)

warm_morph_weight_loss_warm = inner_join(low_food_warm, 
                                         low_food_warm_final, 
                                         by = c('Pop', 
                                                'Rep', 
                                                'Temp')) %>% 
  select(-'low_food_amount') %>% 
  mutate(Weight_dif = avg_weight_final - avg_weight, 
         Length_growth = avg_length_final - avg_length, 
         cond_fac_diff = low_cond_fac_final - low_cond_fac)%>% 
  select(Pop, 
         Rep, 
         Temp,
         Weight_dif, 
         Length_growth, 
         cond_fac_diff) %>% 
  ungroup() %>% 
  filter(Pop %in% c('AWAW', 
                    'MWMW')) %>% 
  summarize(mean_weight_dif = mean(Weight_dif), 
            mean_length_dif = mean(Length_growth), 
            mean_cond_dif = mean(cond_fac_diff))%>% 
  mutate(percent_weight_loss = mean_weight_dif*100, 
         percent_length_change = mean_length_dif*100)

inner_join(low_food_warm, 
           low_food_warm_final, 
           by = c('Pop', 
                  'Rep', 
                  'Temp')) %>% 
  select(-'low_food_amount') %>% 
  mutate(Weight_dif = avg_weight_final - avg_weight, 
         Length_growth = avg_length_final - avg_length, 
         cond_fac_diff = low_cond_fac_final - low_cond_fac)%>% 
  select(Pop, 
         Rep, 
         Temp,
         Weight_dif, 
         Length_growth, 
         cond_fac_diff) %>% 
  ggplot() +
  geom_point(aes(x = Weight_dif, 
                 y = Length_growth, 
                 col = Pop))

# Cold side ---------------------------------------------------------------


low_food_cold = treatment_data %>% 
  # filter(Treatment == 'Low Food', 
  #        Population == 'ACAC') %>% 
  filter(Treatment == 'Low Food', 
         Temperature == '12') %>% 
  group_by(Rep, 
           Population,
           Temperature)%>% 
  rename(Pop = Population,
         Temp = Temperature) %>% 
  summarize(avg_weight = mean(Weight), 
            avg_length = mean(Length)) %>% 
  mutate(low_food_amount = avg_weight*0.02, 
         low_cond_fac = avg_weight/avg_length^(1/3)*100)

# treatment_data %>% 
#   filter(Temperature == '12') %>% 
#   # group_by(Crosses) %>% 
#   distinct(Crosses, 
#            .keep_all = T) %>% 
#   select(Population, 
#          Crosses) %>% 
#   write_csv('Cold_side_Experimental_Crosses.csv')
# treatment_data %>% 
#   filter(Treatment == 'Low Food', 
#          Temperature == '12') %>% 
#   group_by(Population, 
#            Rep) %>% 
#   summarize(Number = n())


high_food_cold = treatment_data %>% 
  # filter(Treatment == 'High Food', 
  #        Population == 'ACAC') %>% 
  filter(Treatment == 'High Food', 
         Temperature == '12') %>% 
  group_by(Rep, 
           Population,
           Temperature)%>% 
  rename(Pop = Population, 
         Temp = Temperature) %>% 
  summarize(avg_weight = mean(Weight), 
            avg_length = mean(Length)) %>% 
  mutate(high_food_amount = avg_weight*0.3, 
         low_cond_fac = avg_weight/avg_length^(1/3)*100)

# treatment_data %>% 
#   filter(Treatment == 'High Food', 
#          Temperature == '12') %>% 
#   group_by(Population, 
#            Rep) %>% 
#   summarize(Number = n())



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
            avg_length_final = mean(Length),
            cond_fac_final = avg_weight_final/avg_length_final^(1/3)*100)

low_food_cold_final = Final_data %>%
  # filter(Treatment == 'Low Food',
  #        Population == 'ACAC') %>%
  filter(Treatment == 'Low Food',
         Temp == '12') %>%
  group_by(Pop,
           Rep, 
           Temp)%>%
  summarize(avg_weight_final = mean(Weight),
            avg_length_final = mean(Length),
            cond_fac_final = avg_weight_final/avg_length_final^(1/3)*100)



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



# cold expectations -------------------------------------------------------
inner_join(low_food_cold, 
           low_food_cold_final, 
           by = c('Pop', 
                  'Rep', 
                  'Temp')) %>% 
  select(-'low_food_amount') %>% 
  mutate(Weight_dif = avg_weight_final - avg_weight, 
         Length_growth = avg_length_final - avg_length, 
         cond_fac_diff = cond_fac_final - low_cond_fac)%>% 
  select(Pop, 
         Rep, 
         Temp,
         Weight_dif, 
         Length_growth, 
         cond_fac_diff) %>% 
  ungroup() %>%
  summarize(mean_weight_dif = mean(Weight_dif), 
            mean_length_dif = mean(Length_growth), 
            mean_cond_dif = mean(cond_fac_diff)) %>% 
  mutate(percent_weight_loss = mean_weight_dif*100, 
         percent_length_change = mean_length_dif*100)


cold_morph_weight_loss_cold = inner_join(low_food_cold, 
                                         low_food_cold_final, 
                                         by = c('Pop', 
                                                'Rep', 
                                                'Temp')) %>% 
  select(-'low_food_amount') %>% 
  mutate(Weight_dif = avg_weight_final - avg_weight, 
         Length_growth = avg_length_final - avg_length, 
         cond_fac_diff = cond_fac_final - low_cond_fac)%>% 
  select(Pop, 
         Rep, 
         Temp,
         Weight_dif, 
         Length_growth, 
         cond_fac_diff) %>% 
  ungroup() %>% 
  filter(Pop %in% c('ACAC', 
                    'MCMC')) %>% 
  summarize(mean_weight_dif = mean(Weight_dif), 
            mean_length_dif = mean(Length_growth), 
            mean_cond_dif = mean(cond_fac_diff)) %>% 
  mutate(percent_weight_loss = mean_weight_dif*100, 
         percent_length_change = mean_length_dif*100)

warm_morph_weight_loss_cold = inner_join(low_food_cold, 
                                         low_food_cold_final, 
                                         by = c('Pop', 
                                                'Rep', 
                                                'Temp')) %>% 
  select(-'low_food_amount') %>% 
  mutate(Weight_dif = avg_weight_final - avg_weight, 
         Length_growth = avg_length_final - avg_length, 
         cond_fac_diff = cond_fac_final - low_cond_fac)%>% 
  select(Pop, 
         Rep, 
         Temp,
         Weight_dif, 
         Length_growth, 
         cond_fac_diff) %>% 
  ungroup() %>% 
  filter(Pop %in% c('AWAW', 
                    'MWMW')) %>% 
  summarize(mean_weight_dif = mean(Weight_dif), 
            mean_length_dif = mean(Length_growth), 
            mean_cond_dif = mean(cond_fac_diff)) %>% 
  mutate(percent_weight_loss = mean_weight_dif*100, 
         percent_length_change = mean_length_dif*100)

inner_join(low_food_cold, 
           low_food_cold_final, 
           by = c('Pop', 
                  'Rep', 
                  'Temp')) %>% 
  select(-'low_food_amount') %>% 
  mutate(Weight_dif = avg_weight_final - avg_weight, 
         Length_growth = avg_length_final - avg_length, 
         cond_fac_diff = cond_fac_final - low_cond_fac)%>% 
  select(Pop, 
         Rep, 
         Temp,
         Weight_dif, 
         Length_growth, 
         cond_fac_diff) %>% 
  ggplot() +
  geom_point(aes(x = Weight_dif, 
                 y = Length_growth, 
                 col = Pop))

##
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
  mutate(high_weight_diff = weight_dif_warm - weight_dif_cold, 
         high_length_diff = Length_dif_warm - Length_dif_cold) %>% 
  dplyr::select(Pop, 
                Rep,
                high_weight_diff, 
                high_length_diff)





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
  mutate(low_weight_diff = weight_dif_warm - weight_dif_cold, 
         low_length_diff = Length_dif_warm - Length_dif_cold) %>% 
  dplyr::select(Pop, 
                Rep,
                low_weight_diff, 
                low_length_diff)


high_compare_results

low_compare_results





# Warm side Second breeding period ----------------------------------------

# Warm side Initial---------------------------------------------------------------
setwd('~/Parsons_Postdoc/Experiment1/')

# setwd('~/Parsons_Postdoc/Experiment1/')

library(tidyverse)

treatment_data_2nd = read_csv('Experiment1_tagging_F1_Breed2.csv')


## quantify condition factor based on that new paper we just read
low_food_warm_2nd = treatment_data_2nd %>% 
  # filter(Treatment == 'Low Food', 
  #        Population == 'ACAC') %>% 
  filter(Treatment == 'Low_food', 
         Temperature == '18') %>% 
  group_by(Population, 
           Rep, 
           Temperature)%>% 
  rename(Pop = Population, 
         Temp = Temperature) %>% 
  summarize(avg_weight = mean(Weight), 
            avg_length = mean(Length), 
            density = n()) %>% 
  mutate(low_food_amount = avg_weight*0.02, 
         low_food_amoung_density = avg_weight*0.02/density,
         low_cond_fac = avg_weight/avg_length^(1/3)*100)

# View(low_food_warm)
# 
# treatment_data_2nd %>% 
#   filter(Treatment == 'Low_food', 
#          Temperature == '18') %>% 
#   group_by(Population, 
#            Rep) %>% 
#   summarize(Number = n())


high_food_warm_2nd = treatment_data_2nd %>% 
  # filter(Treatment == 'High Food', 
  #        Population == 'ACAC') %>% 
  filter(Treatment == 'High_food', 
         Temperature == '18') %>% 
  group_by(Population, 
           Rep, 
           Temperature)%>% 
  rename(Pop = Population, 
         Temp = Temperature) %>% 
  summarize(avg_weight = mean(Weight), 
            avg_length = mean(Length), 
            density = n()) %>% 
  mutate(high_food_amount = avg_weight*0.3,
         high_food_amount_density = (avg_weight*0.3)/density,
         high_cond_fac = avg_weight/avg_length^(1/3)*100)

# treatment_data_2nd %>% 
#   filter(Treatment == 'High_food', 
#          Temperature == '18') %>% 
#   group_by(Population, 
#            Rep) %>% 
#   summarize(Number = n())




# Warm side data analysis -------------------------------------------------

initial_data = read_csv('Experiment1_setup_data.csv') %>% 
  na.omit() 
initial_data$Temperature = as.factor(initial_data$Temperature)
initial_data$Rep = as.factor(initial_data$Rep)

warm_initial = initial_data %>% 
  filter(Temperature == '18')

nest_df = warm_initial %>% 
  group_by(Population,
           Treatment) %>% 
  nest()

rd_model = function(df) {
  lm(Weight ~ Rep, 
     data = df, 
     type = 'II')
}

weight_aov_df = nest_df %>% 
  mutate(mod = map(data, rd_model))

weight_aov_df = weight_aov_df %>% 
  mutate(tidy = map(mod, broom::tidy),
       glance = map(mod, broom::glance), 
       augment = map(mod, broom::augment), 
       rsq = glance %>% map_dbl("adj.r.squared"), 
       pval = glance %>% map_dbl("p.value"))

initial_data %>% 
  filter(Temperature == '18') %>% 
  ggplot(aes(x = Rep, 
             y = Weight))+
  geom_violin(aes(fill = Rep))+
  geom_boxplot(aes(fill = Rep))+
  # facet_grid(~Population)+
  facet_grid(~Treatment + Population)



final_data = read_csv('Experiment1_Sampling_Good.csv') %>% 
  na.omit() 
final_data$Temp = as.factor(final_data$Temp)
final_data$Rep = as.factor(final_data$Rep)



