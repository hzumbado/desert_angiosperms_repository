# Patterns of angiosperm richness
# Section 02
# JSDMs
# Script 02 data preparation

# setup -------------------------------------------------------------------

rm(list = ls())

library(terra)
library(tidyverse)

# rasters -----------------------------------------------------------------

envs <-
  rast('rasters/env_data.tif')

# data --------------------------------------------------------------------

read_rds('data/raw/jsdm/cv_jsdm_data.rds') %>%
  list2env(.GlobalEnv)

# annuals 114 plots, 434 species
# perennials 114 plots, 735 species

annual_data <-
  envs %>%
  terra::extract(
    annual %>%
      dplyr::select(x, y)) %>%
  as_tibble() %>%
  bind_cols(annual) %>%
  select(!ID) %>%
  relocate(plot:y_utm, .before = BIO_03) %>%
  na.omit()

perennial_data <-
  envs %>%
  terra::extract(
    perennial %>%
      dplyr::select(x, y)) %>%
  as_tibble() %>%
  bind_cols(perennial) %>%
  select(!ID) %>%
  relocate(plot:y_utm, .before = BIO_03) %>%
  na.omit()

# two plots removed (NA values)

# save data ---------------------------------------------------------------

list(
  annual_data = annual_data,
  perennial_data = perennial_data) %>%
  write_rds('data/processed/jsdm/jsdm_model_data.rds')
