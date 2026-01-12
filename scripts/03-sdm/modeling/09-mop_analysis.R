# Patterns of angiosperm richness
# Section 03b modeling
# SDMs
# Script 09 MOP analyses

# setup -------------------------------------------------------------------

rm(list = ls())

library(terra)
library(sf)
library(smop)
library(tmap)
library(tidyverse)

# shapefiles --------------------------------------------------------------

list.files(
  'shapefiles',
  pattern = '^cv_w|sa|ca',
  full.names = TRUE) %>%
  map(
    ~ .x %>%
      read_sf() %>%
      st_make_valid) %>%
  set_names(
    'm',
    'cv',
    'salton') %>%
  list2env(.GlobalEnv)

# rasters -----------------------------------------------------------------

# Environmental predictors (present)

present <-
  rast('rasters/env_data.tif') %>%
  crop(m, mask = TRUE)

future_files <-
  list.files(
  'rasters',
  pattern = '\\.tif$',
  full.names = TRUE) %>%
  discard(
    ~ str_detect(.x, '(wc|env)')) %>%
  map(
    ~.x %>%
      rast() %>%
      crop(m, mask = TRUE)) %>%
  set_names(
    'future_126',
    'future_245',
    'future_585')

# Ensure identical variable order

names(present) <- names(future_files$future_126)

stopifnot(
  names(present) == names(future_files$future_126))

# mop analyses ------------------------------------------------------------

mop_ssps <-
  future_files %>%
  map(
    ~ mop(
      M_calibra = present,
      G_transfer = .x,
      percent = 10,
      normalized = TRUE,
      standardize_vars = TRUE))

plot(
  mop_ssps$future_126,
  main = "MOP – SSP126")

extrap_summary <-
  mop_ssps %>%
  imap_dfr(~ {
    vals <- values(.x, mat = FALSE)
    tibble(
      ssp = paste0("SSP", .y),
      total_pixels =
        sum(!is.na(vals)),
      strict_extrap_pixels =
        sum(vals == 0, na.rm = TRUE),
      percent_strict_extrap =
        100 * strict_extrap_pixels / total_pixels)})

strict_map_126 <- mop_ssps$future_126 == 0
strict_map_245 <- mop_ssps$future_245 == 0
strict_map_585 <- mop_ssps$future_585 == 0

plot(
  strict_map_126,
  main = "Strict extrapolation – SSP126")

# map ---------------------------------------------------------------------

map <-
  tm_shape(salton) +
  tm_polygons('darkcyan') +
  tm_shape(m) +
  tm_borders() +
  tm_shape(cv) +
  tm_borders(
    'white',
    lwd = 4)

ssp_126 <-
  strict_map_126 %>%
  tm_shape() +
  tm_raster(
    col.scale =
      tm_scale_categorical(
        values = c('gray70', 'black')),
    col.legend =
      tm_legend(
        title = 'Strict extrapolation',
        position = tm_pos_in('right', 'top'))) +
  map

ssp_245 <-
  strict_map_245 %>%
  tm_shape() +
  tm_raster(
    col.scale =
      tm_scale_categorical(values = c('gray70', 'black')),
    col.legend =
      tm_legend(show = FALSE)) +
  map

ssp_585 <-
  strict_map_585 %>%
  tm_shape() +
  tm_raster(
    col.scale =
      tm_scale_categorical(values = c('gray70', 'black')),
    col.legend =
      tm_legend(show = FALSE)) +
  map

extrap <-
  tmap_arrange(
  ssp_126, ssp_245, ssp_585,
  nrow = 1)

tmap_save(
  extrap,
  'output/figures/fig.s10.jpg',
  dpi = 300)

# save rasters ------------------------------------------------------------

strict_map_126 %>%
  writeRaster('output/models/sdm/predictions/mop/ssp126.tif')

strict_map_245 %>%
  writeRaster('output/models/sdm/predictions/mop/ssp245.tif')

strict_map_585 %>%
  writeRaster('output/models/sdm/predictions/mop/ssp585.tif')
