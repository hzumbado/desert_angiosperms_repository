# Patterns of angiosperm richness
# Section 03a preprocessing
# SDMs
# Script 05 binary polygons list

# Setup -------------------------------------------------------------------

rm(list = ls())

library(sf)
library(tidyverse)

# folders -----------------------------------------------------------------

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

bin_pol %>%
  write_rds(
    paste0(
      polygons,
      'polygons.rds'))
