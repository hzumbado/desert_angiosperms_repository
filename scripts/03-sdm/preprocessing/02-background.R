# Patterns of angiosperm richness
# Section 03a preprocessing
# SDMs
# Script 02 background data

# setup -------------------------------------------------------------------

rm(list = ls())

library(sf)
library(terra)
library(tmap)
library(ecospat)
library(tidyverse)

# raster ------------------------------------------------------------------

envs_stack <-
  rast('rasters/env_data.tif')

# data --------------------------------------------------------------------

occs  <-
  read_csv('data/processed/sdm/model_occs.csv') %>%
  mutate(
    species =
      species %>%
      str_replace(' ', '_'))

my_species <-
  occs %>%
  arrange(species) %>%
  distinct(species) %>%
  pull() %>%
  str_replace(' ', '_')

# sf object ---------------------------------------------------------------

occs_sf <-
  occs %>%
  group_split(species) %>%
  map(
    ~ .x %>%
      st_as_sf(
        coords = c('x', 'y'),
        crs = 4326,
        remove = F)) %>%
  set_names(my_species)

p <-
  occs_sf %>%
  map(
    ~ .x %>%
      vect())

# background --------------------------------------------------------------

model_species <-
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
  # 'Hesperoyucca_whipplei'
  #'Krameria_bicolor'
  #'Langloisia_setosissima'
  #'Larrea_tridentata'
  #'Logfia_depressa'
  #'Nemacladus_glanduliferus'
  #'Oligomeris_linifolia'
  #'Penstemon_centranthifolius'
  #'Psorothamnus_schottii'

background <-
  predicts::backgroundSample(
    mask = envs_stack,
    n = 10000,
    p = p[[model_species]],
    excludep = TRUE) %>%
  as_tibble() %>%
  mutate(
    species = model_species,
    presence = 0) %>%
  bind_rows(
    occs %>%
      filter(species == model_species)) %>%
  write_rds(
    paste0(
    'data/processed/premodel/',
    model_species,
    '_premodel_data.rds'))
