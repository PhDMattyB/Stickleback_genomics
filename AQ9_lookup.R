##############################
## AQ9 stickleback look up
##
## Matt Brachmann (PhDMattyB)
##
## 01.11.2023
##
##############################


library(tidyverse)

setwd('~/Parsons_Postdoc/')

df = read_csv('AQ9_Tank_map_R.csv')

tank_lookup = function(data, Tank_num, num){
  Tank_num = rlang::parse_expr(quo_name(enquo(Tank_num)))
  
  data %>% 
    filter(Tank_num == num) 
  
}

tank_lookup(data = df, 
            T,
            num = '11')


cross_lookup = function(data, Cross_num, num){
  
  Cross_num = rlang::parse_expr(quo_name(enquo(Cross_num)))
  
  data %>% 
    filter(Cross_num == num) 
  
}

cross_lookup(data = df, 
             T,
             num = '64')

df %>% filter(Cross_num == '64')
