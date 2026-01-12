# Patterns of angiosperm richness
# Section 03b modeling
# SDMs
# Script 08 Histogram

# setup -------------------------------------------------------------------

rm(list = ls())

library(terra)
library(sf)
library(tidyverse)

# folders -----------------------------------------------------------------

predictions <- 'output/models/sdm/predictions/average/'

# rasters -----------------------------------------------------------------

present <-
  rast(
    paste0(
      predictions,
    'stacked_x10_present.tif'))

future <-
  rast(
    paste0(
      predictions,
      'stacked_x10_585.tif'))

plot(present)
plot(future)

# plot --------------------------------------------------------------------

x10_present <-
  values(present) %>%
  as_tibble() %>%
  mutate(
    Suitability =
      replace_na(Suitability, 0),
    period = 'present')

x10_future <-
  values(future) %>%
  as_tibble() %>%
  mutate(
    sum =
      replace_na(sum, 0),
    period = 'future') %>%
  rename(Suitability = sum)

histogram <-
  bind_rows(
    x10_present,
    x10_future) %>%
  ggplot(aes(x = Suitability, fill = period)) +
  geom_density(alpha = 0.7) +
  scale_fill_manual(
    values = c('salmon', 'cornflowerblue'),
    labels = c('Future (2081-2100)', 'Present')) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(
    fill = 'Period',
    x = 'Species',
    y = 'Density') +
  theme_classic() +
  theme(legend.position = c(0.8, 0.6))

ggsave('output/figures/histogram.jpg')
