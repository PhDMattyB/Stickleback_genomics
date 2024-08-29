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

# Plot it -----------------------------------------------------------------

ggplot(lipid_clean)+
  geom_violin(aes(x = Season, 
                  y = Fat_content_g, 
                  fill = Season)) +
  geom_point(data = mean_season, 
             aes(x = Season, 
                 y = mean_fat_season), 
             col = 'black', 
             size = 4)

season_temp_col = c('#0077b6',
                   '#e76f51')
ggplot(data = lipid_sum, 
       aes(x = Season, 
           y = mean_fat, 
           group = Lake_morph))+
  ylim(min = 0, 
       max = 0.08)+
  labs(x = 'Season', 
       y = 'Mean fat content')+
  geom_line(col = 'black', 
            size = 1)+
  # geom_pointrange(aes(ymin = mean_fat-sd_fat,
  #                     ymax = mean_fat+sd_fat))+
  # geom_errorbarh(height=.2, 
  #                aes(ymin = mean_fat-sd_fat, 
  #                    ymax = mean_fat-sd_fat))+
  geom_point(aes(col = Temperature), 
             size = 4)+
  scale_color_manual(values = season_temp_col)+
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12), 
        legend.position = 'none')+
  scale_x_discrete(limits = rev)

                   