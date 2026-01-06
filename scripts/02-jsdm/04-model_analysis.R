# Patterns of angiosperm richness
# Section 02
# JSDMs
# Script 04 model fit

# setup -------------------------------------------------------------------

rm(list = ls())

library(Hmsc)
library(tidyverse)

# data --------------------------------------------------------------------

read_rds('data/processed/jsdm/jsdm_model_data.rds') %>%
  list2env(.GlobalEnv)

# switch between annuals and perennials to get two models
# Here model for annuals

Y <-
  #read_rds('output/models/jsdm/annual_species_data.rds') %>%
  read_rds('output/models/jsdm/perennial_species_data.rds') %>%
  pluck('Y')

# load models -------------------------------------------------------------

nChains <- 4
samples <- 1000
thin <- 10
transient <- round(0.5*samples*thin)

filename <-
  paste0(
    #'output/models/jsdm/model_annual_chains_',
    'output/models/jsdm/model_perennial_chains_',
    as.character(nChains),
    '_samples_',
    as.character(samples),
    '_thin_',
    as.character(thin),
    '_transient_',
    as.character(transient), '.rds')

# MCMC convergence --------------------------------------------------------

models <- #change for 2 and 3 respectively
  read_rds(filename)

# # change for 2 and 3 for m2 and m3 respectively

mpost <-
  convertToCodaObject(
    models[[1]]) # model diagnostics/ convergence

preds <-
  computePredictedValues(
    models[[1]]) # model performance

# model fit

MF <-
  evaluateModelFit(
    hM = models[[1]],
    predY = preds) # getting some r2, etc

m1 <- # change for m2 and m3 respectively
  tibble(
    Model = 'mFull', # m1
    #Model = 'mENV', # m2
    #Model = 'mSPACE', # m3

  R2 = MF$TjurR2 %>%
  mean(na.rm = T),

  AUC = MF$AUC %>%
    mean(na.rm = T),

  RMSE = MF$RMSE %>%
    mean(na.rm = T))

model_settings <-
  m1 %>%
  bind_rows(m2) %>%
  bind_rows(m3)

model_settings %>%
  write_csv('output/tables/jsdm/annual_model_fit.csv')
  #write_csv('output/tables/jsdm/perennial_model_fit.csv')

# species metrics ---------------------------------------------------------

# best model: m1

mf_df <-
  tibble(
    species = colnames(models[[1]]$Y),
    R2 = MF$TjurR2,
    AUC = MF$AUC,
    RMSE = MF$RMSE) %>%
    arrange(desc(AUC))

# save model species fit --------------------------------------------------

mf_df %>%
  write_csv('output/tables/jsdm/annual_model_results.csv')
  #write_csv('output/tables/jsdm/perennial_model_results.csv')

# representative species --------------------------------------------------

# annual representative species (best AUC)

annual_10 <-
  mf_df %>%
  slice_max(AUC, n = 10) %>%
  pull(species) %>%
  write_rds('data/processed/jsdm/annual_10.rds')

perennial_10 <-
  mf_df %>%
  slice_max(AUC, n = 10) %>%
  pull(species) %>%
#   write_rds('data/processed/jsdm/perennial_10.rds')

# sdm species ----------------------------------------------------------------

# will be used in section 3 SDMs

mf_df %>%
  slice_max(AUC, n = 10) %>%
  pull(species) %>%
  str_replace_all('_', ' ') %>%
  write_rds('data/processed/sdm/annuals_sdm.rds')
  #write_rds('data/processed/sdm/perennials_sdm.rds')
