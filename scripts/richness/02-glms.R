# Patterns of angiosperm richness
# Section 01
# Diversity patterns of native plants
# Script 02

# setup -------------------------------------------------------------------

rm(list = ls())

library(emmeans)
library(tidyverse)

# data --------------------------------------------------------------------

grid_counts <-  
  read_csv('data/processed/glm_data.csv')

# glms --------------------------------------------------------------------

# model 1

glm1 <- 
  glm(
    richness ~ elevation, 
    family = poisson, 
    data = grid_counts)

summary(glm1)
anova(glm1)

# model 2

glm2 <- 
  glm(
    richness ~ elevation + region, 
    family = poisson, 
    data = grid_counts)

summary(glm2)
anova(glm2)

# model 3

glm3 <- 
  glm(
    richness ~ region, 
    family = poisson, 
    data = grid_counts)

summary(glm3)

anova(glm3)

# best model. Compare AIC

AIC(glm1, glm2, glm3) # glm2 lowest AIC

# abundance

model <- 
  aov(abundance ~ region, 
      data = grid_counts)

summary(model)

# contrasts ---------------------------------------------------------------

em <- emmeans(glm2, 'region')

ct <- 
  contrast(em, 'pairwise', adjust = 'Tukey') %>%
  as_tibble() %>%
  select(!df) 

ct %>%
  mutate(
    across(
      estimate:p.value, 
      ~ round(.x, 2))) 