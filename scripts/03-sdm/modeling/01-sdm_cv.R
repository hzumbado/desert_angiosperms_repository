# Patterns of angiosperm richness
# Section 03b modeling
# SDMs
# Script 01 SDM

# setup -------------------------------------------------------------------

rm(list = ls())

library(ENMeval)
library(terra)
library(tidyterra)
library(tidyverse)

# folders -----------------------------------------------------------------

models <- 'output/models/sdm/files/'

# data --------------------------------------------------------------------

data <-
  read_rds('data/processed/sdm/model_data_sdm.rds')

my_species <-
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

# raster ------------------------------------------------------------------

envs <-
  rast('rasters/env_data.tif')

# model data --------------------------------------------------------------

occs <-
  data %>%
  filter(
    species == my_species,
    presence == 1) %>%
  select(x, y)

bg <-
  data %>%
  filter(
    species == my_species,
    presence == 0) %>%
  select(x, y)

# model -------------------------------------------------------------------

sdm <-
  ENMevaluate(
    occs = occs,
    envs = envs,
    bg = bg,
    tune.args =
      list(fc =
             c('L', 'Q', 'H', 'LQ', 'LH', 'QH', 'LQH'),
           rm = 1:4),
    partitions = 'block',
    algorithm = 'maxent.jar',
    doClamp = TRUE,
    overlap = FALSE,
    taxon.name = my_species,
    parallel = TRUE)

sdm %>%
  write_rds(
    paste0(
      models,
      my_species,
      '_sdm.rds'))

sdm <-
  read_rds(
    paste0(
      models,
      my_species,
      '_sdm.rds'))

# model results -----------------------------------------------------------

model_results <-
  sdm@results %>%
  as_tibble()

# model selection ---------------------------------------------------------

opt.seq <-
  sdm@results %>%
  filter(auc.val.avg == max(auc.val.avg)) %>%
  filter(or.10p.avg == min(or.10p.avg)) %>%
  select(
    tune.args,
    auc = 'auc.train',
    AUC = 'auc.val.avg') #data best model #

opt.seq

# best model --------------------------------------------------------------

bm <-
  sdm@models %>%
  pluck(opt.seq$tune.args[1])

bm_results <-
  bm@results %>%
  as.data.frame() %>%
  rownames_to_column(
    var = 'Parameter') %>%
  as_tibble() %>%
  rename(Value =  V1)

# var contribution --------------------------------------------------------

var_contrib <-
  sdm@variable.importance[[opt.seq$tune.args[1]]] %>%
  as_tibble() %>%
  arrange(desc(percent.contribution))

# response curves --------------------------------------------------------

predicts::partialResponse(
  sdm@models[[opt.seq$tune.args]])

# predictions ------------------------------------------------------------

data <-
  data %>%
  filter(species == my_species)

predictions <-
  data %>%
  mutate(
    prediction =
      as.vector(
        terra::predict(
          bm,
          data,
          type = 'cloglog'))) %>%
  select(
    species:y,
    prediction,
    everything())

# best model settings -----------------------------------------------------

# sdm@results %>%
#   as_tibble() %>%
#   mutate(rm = as_factor(rm)) %>%
#   ggplot(aes(
#     x = fc,
#     y = auc.val.avg,
#     group = rm,
#     col = rm)) +
#   geom_point() +
#   scale_color_manual(values = c(
#     '#999999',
#     '#E69F00',
#     '#56B4E9',
#     'salmon')) +
#   geom_line(aes(col = rm)) +
#   theme_classic()
#
# ggsave(
#   paste0(
#   'output/figures/',
#   my_species,
#   '_best_model_settings_plot.jpg'),
#   dpi = 300)

# save results -----------------------------------------------------------

list(
  model_results = model_results,
  opt.seq = opt.seq,
  best_model = bm,
  best_model_results = bm_results,
  predictions = predictions,
  var_contribution = var_contrib) %>%
  write_rds(
    paste0(
      models,
      my_species,
      '_model_results.rds'))
