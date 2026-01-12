# Patterns of angiosperm richness
# Section 03b modeling
# SDMs
# Script 07b future binary predictions

# Setup -------------------------------------------------------------------

rm(list = ls())

library(sf)
library(terra)
library(tidyverse)

# folders -----------------------------------------------------------------

binary1 <- 'output/models/sdm/binary/future/ssp/1_2.6/'
binary2 <- 'output/models/sdm/binary/future/ssp/2_4.5/'
binary3 <- 'output/models/sdm/binary/future/ssp/5_8.5/'
predictions <- 'output/models/sdm/predictions/average/'
polygons <- 'output/models/sdm/binary/present/polygons/'

# species data ------------------------------------------------------------

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

#stacked predictions
# switch between SSPs
# pattern = '126.tif$', path = binary1
# pattern = '245.tif$', path = binary2
# pattern = '585.tif$', path = binary3

# ssp 1-2.6

x10_raster <-
  list.files(
    binary3,
    pattern = '585.tif$',
    full.names = TRUE) %>%
  map(
    ~.x %>%
      rast()) %>%
  set_names(my_species) %>%
  rast() %>%
  sum()

plot(x10_raster)

x10_raster %>%
  writeRaster(
    paste0(
      predictions,
      'stacked_x10_585.tif'),
    overwrite = TRUE)

# stacked_x10_245.tif
# stacked_x10_585.tif

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


# binary polygon ----------------------------------------------------------
bin_pol <-
  list.files(
    polygons,
    pattern = 'gpkg',
    full.names = TRUE) %>%
  map(
    ~ .x %>%
      read_sf() %>%
      st_make_valid()) %>%
  set_names(my_species)

range_map_pol <-
  bin_pol %>%
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

future_preds_x10 %>%
  names(.) %>%
  walk(~ writeRaster(
    future_preds_x10[[.]],
    paste0(
      future,
      .,
      '_',
      species,
      '_x10.tif'),
    overwrite = TRUE))

future_preds_bin <-
  future_preds_x10 %>%
  map(
    ~.x  %>%
      as.polygons() |>
      st_as_sf() %>%
      filter(Suitability == 1))

future_preds_x10 %>%
  names(.) %>%
  walk(~ write_sf(
    future_preds_bin[[.]],
    paste0(
      binary,
      .,
      '_',
      species,
      '.gpkg'),
    overwrite = TRUE))

plot(future_preds_x10$future_585)
