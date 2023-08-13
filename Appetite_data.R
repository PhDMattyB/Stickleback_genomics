##############################
##  Appetite data 
##
## Matt Brachmann (PhDMattyB)
##
## 13.08.2023
##
##############################


setwd('~/Parsons_Postdoc/Icelandic_workshop/')

library(tidyverse)

theme_set(theme_bw())


appetite_data = read_csv('appetiteData.csv') %>% 
  rename(pellets_eaten = trait)

appetite_data = mutate(.data = appetite_data,
                      Location = as.factor(case_when(
                        pair == 'ASH' ~ 'Áshildarholtsvatn',
                        pair == 'MY' ~ 'Mývatn',
                        pair == 'SKR' ~ 'Sauðárkrókur')))

appetite_data = mutate(.data = appetite_data,
                       Type = as.factor(case_when(
                         source_temp == 'w' ~ 'Geothermal',
                         source_temp == 'c' ~ 'Ambient')))
