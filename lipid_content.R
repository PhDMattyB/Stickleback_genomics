##############################
## Fat content stickleback
##
## Matt Brachmann (PhDMattyB)
##
## 29.08.2024
##
##############################

setwd('~/Parsons_Postdoc/Lipid_starvation_data/')

library(tidyverse)

theme_set(theme_bw())

# Clean it ----------------------------------------------------------------


lipid_data = read_csv('Lipid analysis.csv')

lipid_clean = lipid_data %>% 
  dplyr::select(Population...1,
                Population...2,
         Temperature, 
         Season, 
         W1, 
         W2, 
         Fat_content_g, 
         Percent_fat_content) %>% 
  rename(Population = Population...1, 
         Lake_morph = Population...2) %>% 
  # filter(Population != 'STEIN') %>% 
  filter(Population != 'MYV')

lipid_clean %>% 
  filter(Season == 'Winter') %>% 
  group_by(Population, 
           Temperature) %>% 
  summarize(n = n())


# Analyse it --------------------------------------------------------------


lipid_anova = aov(Fat_content_g ~ Population*Temperature*Season, 
                  data = lipid_clean)

summary(lipid_anova)

# Summarize it ------------------------------------------------------------

lipid_sum = lipid_clean %>% 
  group_by(Lake_morph, 
           Temperature, 
           Season) %>% 
  summarise(mean_fat = mean(Fat_content_g),
            sd_fat = sd(Fat_content_g)) 

mean_season = lipid_clean %>% 
  group_by(Season) %>% 
  summarize(mean_fat_season = mean(Fat_content_g))

mean_season = lipid_clean %>% 
  group_by(Season, 
           Lake_morph) %>% 
  summarize(mean_fat_season = mean(Fat_content_g))

# Plot it -----------------------------------------------------------------

season_temp_col = c('#003049', 
                    '#c1121f')
ggplot(data = lipid_sum, 
       aes(x = Season, 
           y = mean_fat, 
           group = Lake_morph))+
  ylim(min = 0, 
       max = 0.1)+
  labs(x = 'Season', 
       y = 'Mean fat content')+
  geom_violin(data = lipid_clean,
              aes(x = Season,
                  y = Fat_content_g,
                  group = Season), 
              fill = 'white', 
              col = 'grey', 
              size = 2) +
  geom_line(col = 'black',
            size = 1,
            position = position_dodge(0.3))+
  geom_point(aes(col = Temperature),
             size = 4,
             position = position_dodge(0.3))+
  # geom_line(col = 'black', 
  #           size = 1)+
  # geom_point(aes(col = Temperature), 
  #            size = 4)+
  scale_color_manual(values = season_temp_col)+
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12), 
        legend.position = 'none')+
  scale_x_discrete(limits = rev)

ggsave('~/Parsons_Postdoc/Written_things/Stickle_genomic_paper/Fat_Phenotype.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm', 
       width = 15, 
       height = 10)
