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

df %>% 
  filter(Side == 'Cold') %>% 
  filter(Tank_num == '249')


tank_lookup = function(filename, Side = '', Tank_num = ''){
  df = read_csv(filename) 
  
  df = df %>% 
    filter(Side == '') %>% 
    filter(Tank_num == '')
  # 
  # if(Side != 'all'){
  #   sub_df = df %>%
  #     filter(Side == Side)
  #   return(sub_df)
  # }
  # if(Tank_num != 'all'){
  #   sub_df2 = sub_df %>%
  #     filter(Tank_num == Tank_num)
  #   return(sub_df2)
  # }
  
return(df)
  
  
}


tank_lookup(filename = 'AQ9_Tank_map_R.csv', 
            Side = 'Cold', 
            Tank_num = '249')



df %>% 
  filter(Side == 'Cold') %>% 
  filter(Cross_num == '203')
