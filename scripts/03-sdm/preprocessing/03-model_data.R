# Patterns of angiosperm richness
# Section 03a preprocessing
# SDMs
# Script 03 model data prep

# setup -------------------------------------------------------------------

rm(list = ls())

library(sf)
library(terra)
library(tmap)
library(ecospat)
library(fuzzySim)
library(tidyverse)

# raster ------------------------------------------------------------------

envs_stack <-
  rast('rasters/env_data.tif')

# data --------------------------------------------------------------------

background_rarefied <-
  list.files(
  'data/processed/sdm/premodel',
  pattern = 'pre',
  full.names = TRUE) %>%
  map(
    ~.x %>%
      read_rds()) %>%
  bind_rows() %>%
  filter(presence == 0) %>%
  as.data.frame() %>%
  ecospat.occ.desaggregation(
    min.dist = res(envs_stack),
    by = 'species')

background_rarefied %>%
  summarise(
    n = n(),
    .by = species)

background <-
  background_rarefied %>%
  as_tibble()

background %>%
  write_csv('data/processed/sdm/background_rarefied.csv')

occs <-
  list.files(
    'data/processed/sdm/premodel',
    pattern = 'pre',
    full.names = TRUE) %>%
  map(
    ~.x %>%
      read_rds() %>%
      filter(presence == 1)) %>%
  bind_rows() %>%
  bind_rows(background)

occs %>%
  write_csv('data/processed/sdm/model_occs_rarefied.csv')

# grid records ------------------------------------------------------------

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
  #'Hesperoyucca_whipplei'
  #'Krameria_bicolor'
  #'Langloisia_setosissima'
  #'Larrea_tridentata'
  #'Logfia_depressa'
  #'Nemacladus_glanduliferus'
  #'Oligomeris_linifolia'
  #'Penstemon_centranthifolius'
  #'Psorothamnus_schottii'

data <-
  gridRecords(
    rst = envs_stack,
    pres.coords = occs %>%
      filter(presence == 1) %>%
      filter(species == model_species) %>%
      select(x, y),
    abs.coords = occs %>%
      filter(presence == 0) %>%
      filter(species == model_species) %>%
      select(x, y),
    na.rm = T) %>%
  as_tibble() %>%
  mutate(species = model_species) %>%
  select(species, everything()) %>%
  na.omit() %>%
  write_csv(
    paste0(
      'data/processed/sdm/model/',
      model_species,
      '_model_data.csv'))

model_data <-
  list.files(
  'data/processed/sdm/model/',
  pattern = '.csv$',
  full.names = TRUE) %>%
  map(
    ~.x %>%
      read_csv()) %>%
  bind_rows() %>%
  write_rds('data/processed/model_data_sdm.rds')

# map ---------------------------------------------------------------------

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

occs_sf <-
  model_data %>%
  filter(presence == 1) %>%
  st_as_sf(
    coords = c('x', 'y'),
    crs = 4326,
    remove = F)

tm_shape(cv_25) +
  tm_borders(lwd = 2) +
  tm_shape(cv) +
  tm_borders(lwd = 2, 'red') +
  tm_shape(salton) +
  tm_polygons('lightblue') +
  tm_shape(occs_sf) +
  tm_symbols(
    size = 0.3,
    fill = 'species',
    tm_scale_categorical(values = 'brewer.set2')) +
  tm_scalebar(
    position = c('left', 'bottom')) +
  tm_layout(legend.outside = TRUE)
