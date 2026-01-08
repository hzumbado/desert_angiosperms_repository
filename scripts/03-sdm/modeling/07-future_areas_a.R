# Patterns of angiosperm richness
# Section 03b modeling
# SDMs
# Script 07a future binary predictions

# Setup -------------------------------------------------------------------

rm(list = ls())

library(sf)
library(terra)
library(tidyverse)

# folders -----------------------------------------------------------------

models <- 'output/models/sdm/files/'
future <- 'output/models/sdm/predictions/future/'

# rasters -----------------------------------------------------------------

layer_names <-
  list.files(
    'rasters',
    pattern = 'future',
    full.names = FALSE) %>%
  str_remove('.tif')

envs_future <-
  list.files(
    'rasters',
    pattern = 'future',
    full.names = TRUE) |>
  map(~ (rast(.x))) |>
  set_names(layer_names)

# species -----------------------------------------------------------------

# run once for each species

species <-
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

# best model results ------------------------------------------------------

best_model <-
  read_rds(
    paste0(
      models,
      species,
      '_model_results.rds')) %>%
  pluck('best_model')

# future preds ------------------------------------------------------------

future_preds <-
  envs_future  %>%
  map(
    ~ .x %>%
      terra::predict(
        best_model,
        args = c("outputformat=logistic"),
        na.rm = T)) %>%
  set_names(layer_names)

names(future_preds[[1]]) <- 'Suitability'
names(future_preds[[2]]) <- 'Suitability'
names(future_preds[[3]]) <- 'Suitability'

plot(future_preds[[1]])

future_preds %>%
  names(.) %>%
  walk(~ writeRaster(
    future_preds[[.]],
    paste0(
      future,
      .,
      '_',
      species,
      ".tif"),
    overwrite = TRUE))
