# Patterns of angiosperm richness
# Section 03b modeling
# SDMs
# Script 07b future binary rasters and polygons

# Setup -------------------------------------------------------------------

rm(list = ls())

library(sf)
library(terra)
library(tidyverse)

# folders -----------------------------------------------------------------

models <- 'output/models/sdm/files/'
future <- 'output/models/sdm/predictions/future/'
binary_r <- 'output/models/sdm/binary/future/clamped_preds/'
binary_p <- 'output/models/sdm/binary/future/polygons/'
polygons <- 'output/models/sdm/binary/present/polygons/'

# rasters -----------------------------------------------------------------

layer_names <-
  list.files(
    'rasters',
    pattern = 'future',
    full.names = FALSE) %>%
  str_remove('.tif')

# species -----------------------------------------------------------------

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

# x10 ---------------------------------------------------------------------

# x10 threshold

x10 <-
  read_rds(
    paste0(
      models,
      species,
      '_model_results.rds')) %>%
  pluck('best_model_results') %>%
  rename_all(., .funs = tolower) %>%
  filter(
    parameter ==
      'X10.percentile.training.presence.Cloglog.threshold') %>%
  dplyr::select(value) %>%
  pull()

# future preds ------------------------------------------------------------

future_preds <-
  list.files(
    future,
    pattern = species,
    full.names = TRUE) %>%
  map(
    ~ .x %>%
      rast()) %>%
  set_names(layer_names)

# binary polygons ---------------------------------------------------------

# see sdm/preprocessing 05

range_map_pol <-
  read_rds(
    paste0(
      polygons,
      'polygons.rds')) %>%
  pluck(species)

future_preds_x10 <-
  future_preds %>%
  map(
    ~ .x %>%
      crop(
        range_map_pol,
        mask = TRUE) %>%
      clamp(
        lower = x10,
        values = FALSE))

# save binary rasters

future_preds_x10 %>%
  names(.) %>%
  walk(~ writeRaster(
    future_preds_x10[[.]],
    paste0(
      binary_r,
      .,
      '_',
      species,
      '_x10.tif'),
    overwrite = TRUE))

plot(future_preds$future_585)

# save binary polygons

future_preds_bin <-
  future_preds_x10 %>%
  map(
    ~.x  %>%
      as.polygons() |>
      st_as_sf() %>%
      filter(Suitability == 1))

future_preds_bin %>%
  names(.) %>%
  walk(~ write_sf(
    future_preds_bin[[.]],
    paste0(
      binary_p,
      .,
      '_',
      species,
      '.gpkg'),
    overwrite = TRUE))

plot(future_preds_x10$future_585)
