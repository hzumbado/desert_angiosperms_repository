# Patterns of angiosperm richness
# Section 03b modeling
# SDMs
# Script 06 Model evaluation

# setup -------------------------------------------------------------------

rm(list = ls())

library(modEvA)
library(ecospat)
library(tidyverse)

# folders -----------------------------------------------------------------

models <- 'output/models/sdm/files/'

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

# data --------------------------------------------------------------------

model_data <-
  list.files(
    models,
    pattern = 'model_results.rds',
    full.names = T) %>%
  map(~ .x %>%
  read_rds() %>%
    pluck('predictions')) %>%
    set_names(my_species)

opt_seq <-
  list.files(
    models,
    pattern = 'model_results.rds',
    full.names = T) %>%
  map(~ .x %>%
        read_rds() %>%
        pluck('opt.seq')) %>%
  set_names(my_species)

model_results <-
  list.files(
    models,
    pattern = 'model_results.rds',
    full.names = T) %>%
  map(~ .x %>%
        read_rds() %>%
        pluck('best_model_results')%>%
        rename_all(., .funs = tolower) %>%
        filter(
          parameter %in% c(
            'X.Training.samples',
            'X10.percentile.training.presence.Cloglog.threshold'))) %>%
  set_names(my_species)

contribution <-
  list.files(
  models,
  pattern = 'model_results.rds',
  full.names = T) %>%
  map(~ .x %>%
        read_rds() %>%
        pluck('var_contribution') %>%
        arrange(variable)) %>%
  set_names(my_species)

category_fn <-
  function(df, name){

  mutate(df,
         Species = name,
         .before = 1)
  }

envs_contribution <-
  imap(
    contribution,
    ~ category_fn(.x, .y)) %>%
  bind_rows() %>%
  mutate(
    across(
      3:4,
      ~ round(.x, 1)),
  permutation =
    paste0('(', permutation.importance, ')')) %>%
  unite(
    'contribution',
    c(percent.contribution, permutation), sep = ' ') %>%
  select(!permutation.importance) %>%
  pivot_wider(
    names_from = variable,
    values_from = contribution)

envs_contribution %>%
  write_csv('output/tables/sdm/var_contribution.csv')


# THRESHOLD-DEPENDENT MEASURES (classification) ####

sp <-
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

# threshMeasures ----------------------------------------------------------

par(mar = c(6, 3, 2, 1))

measures <-
  threshMeasures(
  obs = model_data[[sp]]$presence,
  pred = model_data[[sp]]$prediction,
  thresh = 'maxSSS',
  main = 'MXT',
  measures = 'TSS')

prev <- measures$Prevalence

eval <-
  measures$ThreshMeasures %>%
  as.data.frame() %>%
  rownames_to_column(
    var = 'Parameter') %>%
  as_tibble() %>%
  pivot_wider(
    names_from = Parameter,
    values_from = Value)

# Boyce index -------------------------------------------------------------

data <-
  model_data[[sp]] %>%
  filter((!is.na(prediction)))

boyce_index <-
  ecospat.boyce(
  fit = data$prediction,
  obs = data$prediction[data$presence == 1])

boyce_index$cor

# eval table --------------------------------------------------------------

tibble(
  species = sp,
  prev,
  eval,
  boyce = boyce_index$cor) %>%
  rename_all(., .funs = tolower)
