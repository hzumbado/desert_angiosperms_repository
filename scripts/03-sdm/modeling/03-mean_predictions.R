# Patterns of angiosperm richness
# Section 03b modeling
# SDMs
# Script 03 SDM mean predictions

# setup -------------------------------------------------------------------

rm(list = ls())

library(terra)
library(tidyterra)
library(sf)
library(tidyverse)

# paths -------------------------------------------------------------------

models <- 'output/models/sdm/files/'
predictions <- 'output/models/sdm/predictions/'
present <- 'output/models/sdm/predictions/present/'

# Species data ------------------------------------------------------------

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

#  shapefiles -------------------------------------------------------------

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

# occurrences -------------------------------------------------------------

occs <-
  read_rds('data/processed/sdm/model_data_sdm.rds') %>%
  filter(presence == 1) %>%
  select(species, x, y) %>%
  st_as_sf(
    coords = c(
      x = 'x',
      y = 'y'),
    crs = 4326) %>%
  st_filter(cv_25) %>%
  vect()

# predictions -------------------------------------------------------------

preds <-
  list.files(
    present,
    full.names = TRUE,
    pattern = '.tif') %>%
  map(
    ~ rast(.x)) %>%
  set_names(my_species) %>%
  rast()

range_map <- mean(preds)

names(range_map) <- 'Mean Suitability'

# save raster -------------------------------------------------------------

range_map %>%
  writeRaster(
    paste0(
      predictions,
      'average/plants_mean_predictions.tif'),
    overwrite = TRUE)

# plot raster -------------------------------------------------------------

rm <-
  range_map %>%
  as.polygons()

plot(range_map)

plot(rm, add = T)

points(
  occs,
  pch = 19,
  cex = 0.5,
  border = 'black',
  col = c('red', 'green', 'blue'))
