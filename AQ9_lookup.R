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


df %>% 
  filter(Side == 'Cold') %>% 
  filter(Cross_num == '203')
