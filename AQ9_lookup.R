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

tank_lookup = function(data, Side, Side_room, Tank_num, num){
  
  # Side <- enquo(Side)
  # Tank_num <- enquo(Tank_num)
  Side = rlang::parse_expr(quo_name(enquo(Side)))
  Tank_num = rlang::parse_expr(quo_name(enquo(Tank_num)))
  
  data %>% 
    filter(Side == Side_room, 
           Tank_num == num) 
  
}


tank_lookup(data = df, 
            Side,
            Side_room = 'Cold', 
            Tank_room,
            num = '250')

# 
# 
# df %>% 
#   filter(Side == 'Cold') %>% 
#   filter(Cross_num == '203')
