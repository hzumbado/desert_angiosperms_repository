# Patterns of angiosperm richness
# Section 02
# JSDMs
# Script 01

# setup -------------------------------------------------------------------

rm(list = ls())

library(terra)
library(usdm)
library(tidyverse)

# climatic data -----------------------------------------------------------

wc <-
  rast('rasters/wc_0.5.tif')

# vif analyses ------------------------------------------------------------

vif <-
  vifstep(wc, th = 3)

vif@results %>%
  as_tibble() %>%
  mutate(
    Variables =
      Variables %>%
      fct_recode(
        'BIO 03' = 'BIO_03',
        'BIO 07' = 'BIO_07',
        'BIO 08' = 'BIO_08',
        'BIO 09' = 'BIO_09'))

envs <-
  wc %>%
  subset(vif@results$Variables)

plot(envs[[1]]) # BIO 03

envs %>%
  writeRaster(
    'rasters/env_data.tif',
    overwrite = TRUE)
