# Patterns of angiosperm richness
# Section 03a preprocessing
# SDMs
# Script 01 occurrence thinning

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

my_species <-
  c(
    'Achyronychia cooperi',
    'Agave deserti',
    'Antirrhinum filipes',
    'Atriplex hymenelytra',
    'Camissoniopsis pallida',
    'Caulanthus lasiophyllus',
    'Cistanthe ambigua',
    'Cylindropuntia ganderi',
    'Eremothera boothii',
    'Ericameria linearifolia',
    'Galium stellatum',
    'Hesperoyucca whipplei',
    'Krameria bicolor',
    'Langloisia setosissima',
    'Larrea tridentata',
    'Logfia depressa',
    'Nemacladus glanduliferus',
    'Oligomeris linifolia',
    'Penstemon centranthifolius',
    'Psorothamnus schottii')

my_species <-
  read_rds(
    'data/raw/sdm/stacked_model_species.rds') %>%
  arrange(species) %>%
  distinct(species) %>%
  pull() %>%
  str_replace(' ', '_')

occs <-
  read_rds(
    'data/raw/sdm/stacked_model_species.rds') %>%
  group_split(species) %>%
  set_names(my_species) %>%
  map(~.x %>%
        st_as_sf(
    coords = c('x', 'y'),
    crs = 4326,
    remove = F)) %>%
  bind_rows() %>%
  as_tibble() %>%
  select(x, y, species) %>%
  as.data.frame() %>%
  ecospat.occ.desaggregation(
    min.dist = res(envs_stack),
    by = 'species')

occs %>%
  as_tibble() %>%
  summarise(
    n = n(),
    .by = species) %>%
  write_csv('output/tables/sdm/model_occs.csv')

occs %>%
  as_tibble() %>%
  mutate(presence = 1) %>%
  relocate(
    species,
    .before = x) %>%
  write_csv('data/processed/sdm/model_occs.csv')
