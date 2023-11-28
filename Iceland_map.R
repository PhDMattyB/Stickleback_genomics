##############################
## iceland map data 
##
## Matt Brachmann (PhDMattyB)
##
## 03.08.2023
##
##############################
setwd('~/Parsons_Postdoc/')

library(tidyverse)
library(maps)

map_data = map_data('world') %>% 
  as_tibble() %>% 
  filter(region == 'Iceland') 

theme_set(theme_minimal())  

iceland_map = ggplot(map_data) +
  geom_map(data = map_data, 
           map = map_data, 
           aes(x = long, 
               y = lat, 
               map_id = region), 
           col = 'black',
           size = 1,
           fill = 'white')+
  theme(panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank())

ggsave('Iceland_map_skinny.tiff', 
       plot = iceland_map, 
       dpi = 'retina')
