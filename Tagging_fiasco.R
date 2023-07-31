##############################
## Fish tagging fiasco
##
## Matt Brachmann (PhDMattyB)
##
## 31.07.2023
##
##############################

library(tidyverse)

setwd('~/Parsons_Postdoc/')

fiasco_data = read_csv('Scarred_fish_data.csv')

mean_vals = fiasco_data %>% 
  group_by(Cross) %>% 
  mutate(mean_weight = mean(weight_gm), 
         mean_length = mean(length_mm)) %>%
  summarize(mean_weight = mean(weight_gm), 
            mean_length = mean(length_mm)) 

## The fucking 'problem' cross is bigger than the 'normals'


# Data analysis -----------------------------------------------------------

## lets do some stats
## Anova on body weight per cross
weight_anova = aov(weight_gm ~ Cross, 
                   data = fiasco_data)
summary(weight_anova)
pairwise.t.test(fiasco_data$weight_gm,
                fiasco_data$Cross, 
                p.adjust.method = 'bonf')
## Body weight does not differ per cross


## Anova of length per cross
length_anova = aov(length_mm ~ Cross, 
                   data = fiasco_data)
summary(length_anova)
pairwise.t.test(fiasco_data$length_mm,
                fiasco_data$Cross, 
                p.adjust.method = 'bonf')
## length does not differ per cross

## manova test to see if they differ in a multivariate framework
dependent_vars = cbind(fiasco_data$weight_gm, 
                       fiasco_data$length_mm)

manova = manova(dependent_vars ~ Cross, 
                data = fiasco_data)

summary(manova)
## NOPE the crosses do NOT differ in a multivariate analysis of 
## length and weight per cross


# Graphs ------------------------------------------------------------------


