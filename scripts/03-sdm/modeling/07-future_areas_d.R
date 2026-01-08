# Patterns of angiosperm richness
# Section 03b modeling
# SDMs
# Script 07d Stacked maps and ESH estimation

# setup -------------------------------------------------------------------

rm(list = ls())

library(sf)
library(terra)
library(tmap)
library(tidyverse)

# folders -----------------------------------------------------------------

predictions <- 'output/models/sdm/predictions/average/'
binary <- 'output/models/sdm/binary/future/polygons/'

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

# shapefiles --------------------------------------------------------------

list.files(
  'shapefiles',
  pattern = '^cv|sa',
  full.names = TRUE) %>%
  map(
    ~ .x %>%
      read_sf() %>%
      st_make_valid) %>%
  set_names(
    'cv_25',
    'cv',
    'salton') %>%
  list2env(.GlobalEnv)

# raster ------------------------------------------------------------------

layer_names <-
  list.files(
    predictions,
    pattern = 'x10',
    full.names = FALSE)

x10_raster <-
  list.files(
    predictions,
    pattern = 'x10',
    full.names = TRUE) %>%
  map(
    ~.x %>%
      rast()) %>%
  set_names(layer_names)

# switch between SSPs

map <-
  x10_raster$stacked_x10_585.tif %>%
  tm_shape() +
  tm_raster(
    col.scale = tm_scale_continuous(
      values = '-kovesi.linear_gow_65_90_c35'),
    col.legend = tm_legend(
      title = 'Species',
      orientation = 'landscape',
      width = 60,
      reverse = FALSE)) +
  tm_shape(salton) +
  tm_polygons('lightblue') +
  tm_shape(cv) +
  tm_borders() +
  tm_shape(cv_25) +
  tm_borders()

# map %>%
#   tmap_save('output/figures/stacked_present.png')
#
# map %>%
#   tmap_save('output/figures/stacked_126.jpg')
#
# map %>%
#   tmap_save('output/figures/stacked_245.jpg')

map %>%
  tmap_save('output/figures/stacked_585.jpg')

# area x10 ----------------------------------------------------------------

# switch between SSPs

bin_pol <-
  list.files(
    binary,
    pattern = '_126', # swich here between 1-2.6 and 2-4.5
    full.names = TRUE) %>%
  map(
    ~ .x %>%
      read_sf() %>%
      st_make_valid()) %>%
  set_names(my_species)

esh <-
  bin_pol %>%
  map(
    ~.x %>%
      st_area() %>%
      units::set_units('km^2')) %>%
  as_tibble() %>%
  pivot_longer(
    cols = 1:20,
    names_to = 'species',
    values_to = 'area')

esh %>%
  write_csv('output/tables/sdm/esh_126.csv')

esh %>%
  write_csv('output/tables/sdm/esh_245.csv')

# for 5.85 0 km2 for Penstemon_centranthifolius

bin_pol <-
  list.files(
    binary,
    pattern = '_585', # swich here.
    full.names = TRUE) %>%
  map(
    ~ .x %>%
      read_sf() %>%
      st_make_valid()) %>%
  set_names(my_species)

bin_pol$Penstemon_centranthifolius <- NULL

esh <-
  bin_pol %>%
  map(
    ~.x %>%
      st_area() %>%
      units::set_units('km^2')) %>%
  as_tibble() %>%
  pivot_longer(
    cols = 1:19,
    names_to = 'species',
    values_to = 'area')

esh %>%
  write_csv('output/tables/sdm/esh_585.csv')

# suitability by species number
# switch between SSPs
# all calibration area

a <-
  x10_raster$stacked_x10_585.tif %>%
  as.polygons() %>%
  st_as_sf() %>%
  st_area() %>%
  units::set_units('km^2') %>%
  as_tibble()

# only Coachella Valley

a_cv <-
  x10_raster$stacked_x10_585.tif %>%
  as.polygons() %>%
  crop(vect(cv)) %>%
  st_as_sf() %>%
  st_area() %>%
  units::set_units('km^2') %>%
  as_tibble()

# areas -------------------------------------------------------------------

# as percentage

areapol <-
  \(raster){

    p <-
      raster %>%
      as.polygons() %>%
      st_as_sf() %>%
      st_area()

    p1 <-
      units::set_units(p, 'km^2') %>%
      as_vector()/sum(p/1000000)/10000

    p1 <- round(p1)
    return(p1)
  }

# all calibration area

(x10_areas <- areapol(x10_raster$stacked_x10_585.tif))

# only Coachella Valley

(x10_areas <- areapol(x10_raster$stacked_x10_585.tif %>% crop(vect(cv))))

