# Patterns of angiosperm richness
# Section 03a preprocessing
# SDMs
# Script 04 future environments

# setup -------------------------------------------------------------------

rm(list = ls())

library(sf)
library(terra)
library(tidyterra)
library(tidyverse)

# data --------------------------------------------------------------------

list.files(
  'shapefiles',
  pattern = '^cv|sa',
  full.names = TRUE) %>%
  map(~ .x %>%
        read_sf() %>%
        st_make_valid) %>%
  set_names('cv_25', 'cv', 'salton') %>%
  list2env(.GlobalEnv)

vars <- c('bio3', 'bio7', 'bio8', 'bio9')

layer_names <-
  c(
    'acm_126',
    'acm_245',
    'acm_585',
    'veg_126',
    'veg_245',
    'veg_585',
    'gem_126',
    'gem_245',
    'gem_585')

# here change to the directory where global rasters are saved

list.files(
  'C:/Users/zumba/Documents/rasters/worldclim/',
  pattern = '\\.tif$',
  full.names = TRUE) %>%
  map(~.x %>%
        rast() %>%
        crop(cv_25, mask = TRUE) %>%
        select(3, 7, 8, 9)) %>%
  set_names(layer_names) %>%
  list2env(.GlobalEnv)

names(acm_126) <- vars
names(acm_245) <- vars
names(acm_585) <- vars
names(veg_126) <- vars
names(veg_245) <- vars
names(veg_585) <- vars
names(gem_126) <- vars
names(gem_245) <- vars
names(gem_585) <- vars

future_126 <- mean(acm_126, gem_126, veg_126)
future_245 <- mean(acm_245, gem_245, veg_245)
future_585 <- mean(acm_585, gem_585, veg_585)

varnames(future_126) <- 'mean_126'
varnames(future_245) <- 'mean_245'
varnames(future_585) <- 'mean_585'

plot(future_126$bio3)

envs_cropped <-
  list(
    future_126 = future_126,
    future_245 = future_245,
    future_585 = future_585)

envs_cropped %>%
  names(.) %>%
  walk(~ writeRaster(
    envs_cropped[[.]],
    paste0('rasters/', ., '.tif'),
    overwrite = TRUE))
