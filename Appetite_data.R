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
                         source_temp == 'W' ~ 'Geothermal',
                         source_temp == 'C' ~ 'Ambient')))

appetite_means = appetite_data %>% 
  group_by(Location, 
           Type) %>% 
  summarise(mean_pellets = mean(pellets_eaten))



# statistical test --------------------------------------------------------
## geothermal fish eat more pellets than ambient fish
aov_test = aov(pellets_eaten ~ pair*source_temp, 
               data = appetite_data)

summary(aov_test)


# plot it up --------------------------------------------------------------

appetite_cols = c('#003049', 
                  '#c1121f')

appetite_plot = ggplot(data = appetite_data, 
         aes(x = Location, 
             y = pellets_eaten, 
             fill = Type), 
         col = 'black') +
    # geom_boxplot()+
  geom_violin() +
    stat_summary(fun = "mean", 
                 geom = "point", 
                 position = position_dodge(0.9), 
                 size = 3)+
  # geom_point(data = appetite_means, 
  #            aes(x = Location, 
  #                y = mean_pellets), 
  #            col = 'black', 
  #            size = 3)+
  scale_fill_manual(values = appetite_cols)+
  labs(y = 'Pellets eaten to satiation')+
    theme(panel.grid = element_blank(), 
          axis.title.x = element_blank(), 
          axis.title.y = element_text(size = 14), 
          axis.text = element_text(size = 12))


ggsave('Appetite_difference.tiff', 
       plot = appetite_plot, 
       dpi = 'retina', 
       units = 'cm', 
       height = 10, 
       width = 15)




