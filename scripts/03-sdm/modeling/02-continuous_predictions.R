# Patterns of angiosperm richness
# Section 03b modeling
# SDMs
# Script 02 SDM continuous predictions

# setup -------------------------------------------------------------------

rm(list = ls())

library(terra)
library(sf)
library(tidyterra)
library(tmap)
library(tidyverse)

# paths -------------------------------------------------------------------

processed <- 'data/processed/sdm/'
models <- 'output/models/sdm/files/'
predictions <- 'output/models/sdm/predictions/'

# shapefiles --------------------------------------------------------------

list.files(
  'shapefiles',
  pattern = '^cv|sa',
  full.names = TRUE) %>%
  map(
    ~ .x %>%
      read_sf() %>%
      st_make_valid) %>%
  set_names('cv_25', 'cv', 'salton') %>%
  list2env(.GlobalEnv)

# raster ------------------------------------------------------------------

envs <-
  rast('rasters/env_data.tif')

names(envs) <- c('bio3', 'bio7', 'bio8', 'bio9')

# selecting best model ----------------------------------------------------

my_species <-
  'Achyronychia_cooperi'
  #'Agave_deserti'
  #'Antirrhinum_filipes'
  #'Atriplex_hymenelytra'
  #'Camissoniopsis_pallida'
  #'Caulanthus_lasiophyllus'
  #'Cistanthe_ambigua'
  #'Cylindropuntia_ganderi'
  #'Eremothera_boothii'
  #'Ericameria_linearifolia'
  #'Galium_stellatum'
  #'Hesperoyucca_whipplei'
  #'Krameria_bicolor'
  #'Langloisia_setosissima'
  #'Larrea_tridentata'
  #'Logfia_depressa'
  #'Nemacladus_glanduliferus'
  #'Oligomeris_linifolia'
  #'Penstemon_centranthifolius'
  #'Psorothamnus_schottii'

best_model <-
  read_rds(
    paste0(
      models,
      my_species,
      '_model_results.rds')) %>%
  pluck('best_model')

# raster predictions ------------------------------------------------------

range_map <-
  predicts::predict(
    best_model,
    envs,
    args = c("outputformat=logistic")) %>%
  mask(salton, inverse = TRUE)

names(range_map) <- 'Suitability'

range_map %>%
  writeRaster(
    paste0(
      predictions,
      my_species,
      '_predictions.tif'),
    overwrite = TRUE)

# map ---------------------------------------------------------------------

range_map %>%
  tm_shape() +
  tm_grid(lines = F) +
  tm_raster(
    col.scale =
      tm_scale_continuous(
        values = "brewer.blues",
        limits = c(0, 1)),
    col.legend =
      tm_legend(,
        bg.color = "gray",
        bg.alpha = 0.5,
        reverse = TRUE,
        orientation = 'landscape')) +
  tm_shape(cv) +
  tm_borders(
    col = 'black',
    lwd = 2) +
  tm_shape(salton) +
  tm_polygons('lightblue') +
  tm_shape(
    cv_25) +
  tm_borders() +
  tm_shape()
