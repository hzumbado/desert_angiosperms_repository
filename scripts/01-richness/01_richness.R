# Patterns of angiosperm richness
# Section 01
# Diversity patterns of native plants
# Script 01

# setup -------------------------------------------------------------------

rm(list = ls())

library(sf)
library(tmap)
library(terra)
library(phyloregion)
library(tidyverse)

# shapefiles --------------------------------------------------------------

list.files(
  'shapefiles',
  pattern = '(eco|cv_2|cv_w|salton)',
  full.names = TRUE) %>%
  map(
    ~.x %>%
      read_sf() %>%
      st_transform(crs = 32611) %>%
      st_make_valid) %>%
  set_names(
    'cv_25',
    'cv',
    'ecoregions',
    'salton') %>%
  list2env(.GlobalEnv)

# data --------------------------------------------------------------------

plots_native <-
  read_rds('data/raw/cv_plots_updated.rds') %>%
  pluck('plot_native')

plots_invasive <-
  read_rds('data/raw/cv_plots_updated.rds') %>%
  pluck('plot_invasive')

# 5 more abundant native species

plots_native %>%
  pivot_longer(
    7:ncol(plots_native),
    names_to = 'species',
    values_to = 'value') %>%
  summarise(
    n = sum(value),
    .by = species) %>%
  slice_max(n = 5, n)

# 5 more abundant invasive species

plots_invasive %>%
  pivot_longer(
    7:ncol(plots_invasive),
    names_to = 'species',
    values_to = 'value') %>%
  summarise(
    n = sum(value),
    .by = species) %>%
  slice_max(n = 5, n)

# sparse matrix -----------------------------------------------------------

# native

plot_data <-
  plots_native %>%
  select(7:ncol(plots_native)) %>%
  as.data.frame()

rownames(plot_data) <-
  plots_native %>%
  pull(grid_id)

sparse_matrix <-
  dense2sparse(plot_data)

# invasive

plot_data_invasive <-
  plots_invasive %>%
  select(7:ncol(plots_invasive)) %>%
  as.data.frame()

rownames(plot_data_invasive) <-
  plots_invasive %>%
  pull(grid_id)

# endemism ----------------------------------------------------------------

endemism <- weighted_endemism(sparse_matrix)
coldspots <- coldspots(endemism) # cold spots
hotspots <- hotspots(endemism)  # hotspots

counts_native <-
  tibble(
    grid_id = row.names(plot_data),
    abundance = Matrix::rowSums(plot_data),
    richness = Matrix::rowSums(plot_data > 0),
    coldspots = coldspots,
    hotspots = hotspots)

counts_invasive <-
  tibble(
    grid_id = row.names(plot_data_invasive),
    abundance = Matrix::rowSums(plot_data_invasive),
    richness = Matrix::rowSums(plot_data_invasive > 0))

# grid counts -------------------------------------------------------------

# native ------------------------------------------------------------------

grid_counts_sf <-
  plots_native %>%
  select(grid_id, grid, elevation, region) %>%
  inner_join(counts_native) %>%
  st_as_sf()

# average metrics

grid_counts_sf %>%
  st_drop_geometry() %>%
  mutate(region = as_factor(region)) %>%
  summarise(
    Plots = n(),
    Richness = mean(richness, na.rm = TRUE),
    sd_r = sd(richness, na.rm = TRUE),
    Abundance = mean(abundance, na.rm = TRUE),
    sd_abundance = sd(abundance, na.rm = TRUE))

# average metrics by ecoregion

grid_counts_sf %>%
  st_drop_geometry() %>%
  mutate(region = as_factor(region)) %>%
  summarise(
    Plots = n(),
    Richness = mean(richness, na.rm = TRUE),
    sd_r = sd(richness, na.rm = TRUE),
    Abundance = mean(abundance, na.rm = TRUE),
    sd_abundance = sd(abundance, na.rm = TRUE),
    .by = region)

grid_counts_sf %>%
  st_drop_geometry() %>%
  write_csv('data/processed/glm_data.csv')

# invasives----------------------------------------------------------------

grid_counts_invasive_sf <-
  plots_invasive %>%
  select(grid_id, grid, elevation, region) %>%
  inner_join(counts_invasive) %>%
  st_as_sf()

# average metrics

grid_counts_invasive_sf %>%
  st_drop_geometry() %>%
  mutate(region = as_factor(region)) %>%
  summarise(
    Plots = n(),
    Richness = mean(richness, na.rm = TRUE),
    sd_r = sd(richness, na.rm = TRUE),
    Abundance = mean(abundance, na.rm = TRUE),
    sd_abundance = sd(abundance, na.rm = TRUE))

# average metrics by ecoregion

grid_counts_invasive_sf %>%
  st_drop_geometry() %>%
  mutate(region = as_factor(region)) %>%
  summarise(
    Plots = n(),
    Richness = mean(richness, na.rm = TRUE),
    sd_r = sd(richness, na.rm = TRUE),
    Abundance = mean(abundance, na.rm = TRUE),
    sd_abundance = sd(abundance, na.rm = TRUE),
    .by = region)

# hotspots ----------------------------------------------------------------

richness_plot1 <-
  grid_counts_sf %>%
  rename(Region = 'region') %>%
  ggplot(aes(
    x = Region,
    y = richness,
    col = Region,
    fill = Region)) +
  geom_boxplot(
    show.legend = FALSE,
    col = 'black',
    linewidth = 0.75)  +
  scale_fill_brewer(palette = "Set1") +
  scale_x_discrete(
    labels = c(
      'Colorado desert' = 'Colorado Desert',
      'Mojave desert' = 'Mojave Desert',
      'SoCal mountains' = 'Southern California Mountains',
      'Sonoran desert' = 'Sonoran Desert')) +
  scale_y_continuous(
    expand = c(0,0)) +
  labs(
    x = 'Ecoregion',
    y = 'Richness') +
  theme_classic() +
  theme(
    axis.text.y =
      element_text(size = 12, colour = 'black'),
    axis.text.x =
      element_blank(),,
    axis.title.y = element_text(size = 14),
    axis.title.x = element_blank())

richness_plot2 <-
  grid_counts_sf %>%
  rename(Region = 'region') %>%
  ggplot(aes(
    x = Region,
    y = abundance,
    col = Region,
    fill = Region)) +
  geom_boxplot(
    show.legend = FALSE,
    col = 'black',
    linewidth = 0.75)  +
  scale_fill_brewer(palette = "Set1") +
  scale_y_continuous(
    expand = c(0,0), limits = c(0, 2000)) +
  labs(
    x = 'Ecoregion',
    y = 'Abundance') +
  scale_x_discrete(
    labels = c(
      'Colorado desert' = 'Colorado Desert',
      'Mojave desert' = 'Mojave Desert',
      'SoCal mountains' = 'Southern California Mountains',
      'Sonoran desert' = 'Sonoran Desert')) +
  theme_classic() +
  theme(
    axis.text =
      element_text(size = 12, colour = 'black'),
    axis.title.y = element_text(size = 14),
    axis.title.x = element_blank())



cowplot::plot_grid(
  richness_plot1,
  richness_plot2,
  ncol = 1,
  labels = c('A)', 'B)'),
  scale = 0.9)

ggsave(
  'output/figures/richness_ecoregion.jpg',
  dpi = 300)

# maps richness -----------------------------------------------------------

# native plants

m1 <-
  tm_shape(grid_counts_sf) +
  tm_polygons(
    fill = 'richness',
    fill.scale =
      tm_scale_continuous(
        values = 'brewer.greens'),
    fill_alpha = 0.75,
    fill.legend =
      tm_legend(
        'Richness',
        orientation = 'landscape',
        width = 60,
        title.size = 2)) +
  tm_shape(salton) +
  tm_polygons('lightblue') +
  tm_shape(cv_25) +
  tm_borders('black', lwd = 2) +
  tm_shape(cv) +
  tm_borders('red', lwd = 3) +
  tm_credits(
    'A)',
    position = c('right', 'top'),
    size = 1) +
  tm_layout(
    legend.outside = T)

m2 <-
  tm_shape(grid_counts_sf) +
  tm_polygons(
    fill = 'abundance',
    fill.scale =
      tm_scale_continuous(
        values = 'brewer.blues'),
    fill_alpha = 0.75,
    fill.legend =
      tm_legend(
        'Abundance',
        orientation = 'landscape',
        width = 60,
        title.size = 2)) +
  tm_shape(salton) +
  tm_polygons('lightblue') +
  tm_shape(cv_25) +
  tm_borders('black', lwd = 2) +
  tm_shape(cv) +
  tm_borders('red', lwd = 3) +
  tm_credits(
    'B)',
    position = c('right', 'top'),
    size = 1)
  tm_layout(legend.outside = T)

# invasive plants

m3 <-
  tm_shape(grid_counts_invasive_sf) +
  tm_polygons(
    fill = 'richness',
    fill.scale =
      tm_scale_continuous(
        values = 'brewer.greens'),
    fill_alpha = 0.75,
    fill.legend =
      tm_legend(
        'Richness',
        orientation = 'landscape',
        width = 60,
        title.size = 2)) +
  tm_shape(salton) +
  tm_polygons('lightblue') +
  tm_shape(cv_25) +
  tm_borders('black', lwd = 2) +
  tm_shape(cv) +
  tm_borders('red', lwd = 3) +
  tm_credits(
    'C)',
    position = c('right', 'top'),
    size = 1)
  tm_layout(legend.outside = T)

m4 <-
  tm_shape(grid_counts_invasive_sf) +
  tm_polygons(
    fill = 'abundance',
    fill.scale =
      tm_scale_continuous(
        values = 'brewer.blues'),
    fill_alpha = 0.75,
    fill.legend =
      tm_legend(
        'Abundance',
        orientation = 'landscape',
        width = 60,
        title.size = 2)) +
  tm_shape(salton) +
  tm_polygons('lightblue') +
  tm_shape(cv_25) +
  tm_borders('black', lwd = 2) +
  tm_shape(cv) +
  tm_borders('red', lwd = 3) +
  tm_credits(
    'D)',
    position = c('right', 'top'),
    size = 1)
  tm_layout(legend.outside = T)

map_richness <-
  tmap_arrange(
    m1,  m2, m3, m4,
    ncol = 2)

tmap_save(
  map_richness,
  'output/figures/maps_richness.jpg')


# map hotspots ------------------------------------------------------------

map_hotspots <-
  tm_shape(
    grid_counts_sf,
    is.main = TRUE) +
  tm_borders() +
  tm_shape(
    grid_counts_sf %>%
      filter(
        coldspots == 1)) +
  tm_polygons(fill = 'blue') +
  tm_shape(
    grid_counts_sf %>%
      filter(
        hotspots == 1)) +
  tm_polygons(
    fill = 'red') +
  tm_shape(salton) +
  tm_polygons('lightblue') +
  tm_shape(cv_25) +
  tm_borders('black', lwd = 2) +
  tm_shape(cv) +
  tm_borders('gray50', lwd = 3) +
  tm_add_legend(
    type = 'polygons',
    fill = 'red',
    labels = 'Hotspot',
    position = c('left', 'bottom')) +
  tm_add_legend(
    type = 'polygons',
    fill = 'blue',
    labels = 'Coldspot',
    position = c('left', 'bottom')) +
  tm_layout(legend.text.size = 1)


tmap_save(
  map_hotspots,
  'output/figures/map_hotspots.jpg')
