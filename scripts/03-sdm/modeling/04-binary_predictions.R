# Patterns of angiosperm richness
# Section 03b modeling
# SDMs
# Script 04 Binary predictions

# setup -------------------------------------------------------------------

rm(list = ls())

library(sf)
library(terra)
library(tidyverse)

# folders -----------------------------------------------------------------

predictions <- 'output/models/sdm/predictions/present/'
binary_1 <- 'output/models/sdm/binary/present/polygons/'
binary_2 <- 'output/models/sdm/binary/present/clamped_preds/'
models <- 'output/models/sdm/files/'

# Species data ------------------------------------------------------------

my_species <-
  c(
    'Achyronychia_cooperi',
    'Agave_deserti',
    'Antirrhinum_filipes',
    'Atriplex_hymenelytra',
    'Camissoniopsis_pallida',
    'Caulanthus_lasiophyllus',
    'Cistanthe_ambigua',
    'Cylindropuntia_ganderi',
    'Eremothera_boothii',
    'Ericameria_linearifolia',
    'Galium_stellatum',
    'Hesperoyucca_whipplei',
    'Krameria_bicolor',
    'Langloisia_setosissima',
    'Larrea_tridentata',
    'Logfia_depressa',
    'Nemacladus_glanduliferus',
    'Oligomeris_linifolia',
    'Penstemon_centranthifolius',
    'Psorothamnus_schottii')

# raster ------------------------------------------------------------------

range_map <-
  list.files(
    predictions,
    pattern = '.tif',
    full.names = T) %>%
  map(
    ~.x %>%
      rast()) %>%
  set_names(my_species)

# best model results ------------------------------------------------------

bm_results <-
  list.files(
    models,
    pattern = 'model_results.rds',
    full.names = T) %>%
  map(~ .x %>%
        read_rds() %>%
        pluck('best_model_results') %>%
        rename_all(., .funs = tolower)) %>%
  set_names(my_species)

# x10 ---------------------------------------------------------------------

x10 <-
  bm_results %>%
  map(
    ~ .x %>%
      filter(
        parameter ==
          'X10.percentile.training.presence.Cloglog.threshold') %>%
      dplyr::select(value) %>%
      pull()) %>%
  as_vector()

# suitability maps --------------------------------------------------------

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

suitability <-
  range_map %>%
  pluck(species) %>%
  clamp(
    lower =  x10[species],
    values = TRUE) # values = F

names(suitability) = 'x10'

plot(suitability)

# save suitability maps ---------------------------------------------------

# binary layer ------------------------------------------------------------

#x10

bin_lim <-
  range_map %>%
  pluck(species) >= x10[species]

x10_f <-
  bin_lim %>%
  terra::as.factor()

x10_f %>%
  writeRaster(
    paste0(
      binary_2,
      species,
      '_binary_predictions.tif'),
    overwrite = TRUE)

bin_pol <-
  as.polygons(bin_lim) %>%
  st_as_sf() %>%
  filter(Suitability == 1)

plot(vect(bin_pol))

bin_pol %>%
  write_sf(
    paste0(
      binary_1,
      species,
      '_binary_predictions.gpkg'))
