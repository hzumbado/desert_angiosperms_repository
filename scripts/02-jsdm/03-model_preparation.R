# Patterns of angiosperm richness
# Section 02
# JSDMs
# Script 03 model preparation

# setup -------------------------------------------------------------------

rm(list = ls())

library(tidyverse)
library(Hmsc)

# data --------------------------------------------------------------------

read_rds('data/processed/model_data.rds') %>%
  list2env(.GlobalEnv)

# switch between annuals and perennials to get two models
# Here model for annuals

data <-
  annual_data
  #perennial_data

y <-
  data %>%
  select(11:ncol(data)) %>%
  as.matrix()

XData <-
  data %>%
  mutate(
    region = region %>%
      fct_recode(
        scm = 'SoCal mountains',
        md = 'Mojave desert',
        sd = 'Sonoran desert',
        cd = 'Colorado desert')) %>%
  tibble::column_to_rownames('plot') %>%
  mutate(
    region =
      as.factor(region)) %>%
  as.data.frame()

# prevalence --------------------------------------------------------------

prevalence <-
  colSums(y) %>%
  as_tibble(rownames = 'species') %>%
  rename(occurrences = value) %>%
  mutate(prevalence = occurrences / 112 * 100) %>%
  arrange(desc(prevalence))

# species data ------------------------------------------------------------

# 100 most prevalent species

model_species <-
  prevalence %>%
  slice_head(n = 100) %>%
  pull(species)

Y <- y[, model_species] # subset model species plot data

list(
  prevalence = prevalence,
  model_species = model_species,
  Y = Y) %>%
  write_rds('output/models/jsdm/annual_species_data.rds')
  #write_rds('output/models/perennial_species_data.rds')

Y %>%
  as.data.frame() %>%
  as_tibble() %>%
  pivot_longer(
    1:100,
    names_to = 'species') %>%
  mutate(species = str_replace(species, '_', ' ')) %>%
  arrange(species) %>%
  distinct(species) %>%
  write_csv('output/tables/annuals.csv')
  #write_csv('output/tables/perennials.csv')



# model data --------------------------------------------------------------

XFormula <-
  ~ region +
  BIO_03 +
  BIO_07 +
  BIO_08 +
  BIO_09

studyDesign <-
  data.frame(
    plot =
      as.factor(data$plot))

xy <-
  data %>%
  select(x, y) %>%
  as.matrix()

rownames(xy) <-
  studyDesign[,1]

rL <-
  HmscRandomLevel(sData = xy)

# models ------------------------------------------------------------------

mFull <-
  Hmsc(
    Y = Y,
    XData = XData,
    XFormula = XFormula,
    distr = 'probit',
    studyDesign = studyDesign,
    ranLevels = list('plot' = rL))

mENV <-
  Hmsc(
    Y = Y,
    XData = XData,
    XFormula = XFormula,
    distr = 'probit',
    studyDesign = studyDesign)

mSPACE <-
  Hmsc(
    Y = Y,
    XData = XData,
    XFormula = ~1,
    distr = 'probit',
    studyDesign = studyDesign,
    ranLevels = list('plot' = rL))

head(mFull$X)
head(mENV$X)
head(mSPACE$X)

# model settings ----------------------------------------------------------

nChains <- 4
samples <- 1000
thin <- 10
transient <- round(0.5*samples*thin)

mFull <-
  sampleMcmc(
    mFull,
    thin = thin,
    samples = samples,
    transient = transient,
    nChains = nChains)

mENV <-
  sampleMcmc(
    mENV,
    thin = thin,
    samples = samples,
    transient = transient,
    nChains = nChains)

mSPACE <-
  sampleMcmc(
    mSPACE,
    thin = thin,
    samples = samples,
    transient = transient,
    nChains = nChains)

# save models -------------------------------------------------------------

models <-
  list(mFull, mENV, mSPACE)

filename <-
  paste0(
    'output/models/model_annual_chains_',
    #'output/models/model_perennial_chains_',
    as.character(nChains),
    '_samples_',
    as.character(samples),
    '_thin_',
    as.character(thin),
    '_transient_',
    as.character(transient), '.rds')

filename

models %>%
  write_rds(filename)
