# Patterns of angiosperm richness
# Section 02
# JSDMs
# Script 05 model results

# setup -------------------------------------------------------------------

rm(list = ls())

library(Hmsc)
library(tidyverse)
#library(RColorBrewer)

# data --------------------------------------------------------------------

read_rds('data/processed/model_data.rds') %>%
  list2env(.GlobalEnv)

# switch between annuals and perennials to get two models
# Here model for annuals

Y <-
  read_rds('output/models/annual_species_data.rds') %>%
  #read_rds('output/models/perennial_species_data.rds') %>%
  pluck('Y')

annual_10 <-
  read_rds('data/processed/annual_10.rds')

# perennial_10 <-
#   read_rds('data/processed/perennial_10.rds')

# perennial_10 %>%
#   sort()

# load models -------------------------------------------------------------

nChains <- 4
samples <- 1000
thin <- 10
transient <- round(0.5*samples*thin)

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

# BM convergence --------------------------------------------------------

best_model <-
  read_rds(filename) %>%
  pluck(1)

mpost <-
  convertToCodaObject(
    best_model) # model diagnostics/ convergence

preds <-
  computePredictedValues(
    best_model) # model performance

# model fit ---------------------------------------------------------------

# sample size

ess.beta <-
  effectiveSize(mpost$Beta) %>%
  as_tibble() %>%
  rename(ess_beta = value)

ess.beta$ess_beta %>%
  mean(na.rm = T) # 1357 annual

# beta parameter

psrf.beta <-
  gelman.diag(
    mpost$Beta,
    multivariate = FALSE)$psrf %>%
  as_tibble() %>%
  rename(psrf_beta = `Point est.`)

psrf.beta$psrf_beta %>%
  mean(na.rm = T)

# 1.002 annual
# 1.002 perennial

ess.beta %>%
  bind_cols(psrf.beta)

diag_all2 <-
  ggpubr::ggarrange(
    ggplot(
      ess.beta,
      aes(x = ess_beta)) +
      geom_histogram(
        fill = 'cornflowerblue',
        col = 'black',
        bins = 20) +
      xlab('Effective Sample Size') +
      ylab('Count') +
      scale_x_continuous(expand = c(0,0)) +
      scale_y_continuous(expand = c(0,0)) +
      theme_classic(),
    ggplot(
      psrf.beta, aes(
        x = psrf_beta)) +
      geom_histogram(
        fill = 'cornflowerblue',
        col = 'black',
        bins = 10) +
      xlab('Gelman Diagnostic') +
      ylab('Count') +
      scale_x_continuous(
        expand = c(0,0)) +
      scale_y_continuous(expand = c(0,0)) +
      theme_classic(),
    align = 'v') +
  ggtitle('All Plots')

ggsave(
  diag_all,
  filename = 'output/figures/geldman_ess_m1_annual.jpg',
  #filename = 'output/figures/geldman_ess_m1_perennial.jpg',
  width = 5.5,
  height = 3.5,
  bg = 'white')

# variance partitioning ---------------------------------------------------

VP <-
  computeVariancePartitioning(best_model) # variance partitioning

vp_df <-
  VP$vals %>%
  as_tibble(rownames = 'variable') %>%
  pivot_longer(
    cols = names(.)[2:ncol(.)],
    names_to = 'species',
    values_to = 'value') %>%
  na.omit() %>%
  arrange(species) %>%
  filter(species %in% annual_10) %>%
  #filter(species %in% perennial_10) %>%
  mutate(
    species = species %>%
      str_replace('_', ' '))

# summary

vp_summary <-
  vp_df %>%
  summarise(
    value = mean(value),
    .by = variable) %>%
  arrange(desc(value))

vp_order <-
  vp_df %>%
  filter(
    variable == 'region') %>%
  arrange(value) %>%
  mutate(
    species_f =
      factor(species, levels = .$species)) %>%
  select(species, species_f)

vp <-
  inner_join(
    vp_df,
    vp_order) %>%
  mutate(
    variable = factor(
      variable,
      c( 'region',
         'BIO_03' ,
         'BIO_07',
         'BIO_08'  ,
         'BIO_09',
         'Random: plot'))) %>%
  filter(variable != 'Random: plot' ) %>%
  mutate(
    variable = variable %>%
      fct_recode(
        Region = 'region',
        'BIO 03' = 'BIO_03',
        'BIO 07' = 'BIO_07',
        'BIO 08' = 'BIO_08',
        'BIO 09' = 'BIO_09')) %>%
  ggplot(
    aes(
      x = value,
      y = species_f,
      fill = variable)) +
  geom_bar(stat = 'identity') +
  scale_x_continuous(expand = c(0,0)) +
  theme_classic() +
  ylab('Species') +
  xlab('Proportion of variance explained') +
  scale_fill_brewer(
    name = 'Predictor',
    palette = 'Dark2') +
  theme(
    axis.title =
      element_text(size = 16, colour = 'black'),
    axis.title.x = element_blank(),
    axis.text.y = element_text(face = 'italic', size = 14, colour = 'black'),
    axis.text.x = element_text(size = 14, colour = 'black'),
    legend.position = c(0.9, 0.02),
    legend.text = ggtext::element_markdown(),
    legend.justification = c(0.5,0),
    legend.background = element_rect(color = 'black'),
    legend.key.size = unit(0.5, "cm"), # Adjust the overall key size
    legend.key.width = unit(0.25, "cm"), # Adjust key width
    legend.key.height = unit(0.25, "cm"))
    #legend.position = "none")

vp

ggsave(
  filename = 'output/figures/variance_annual_occurrence.png',
  #filename = 'output/figures/variance_perennial_occurrence.png',
  height = 10,
  width = 15)

# categorical gradient ----------------------------------------------------

gradient <-
  constructGradient(
    best_model,
    focalVariable = 'region')

predY_region <-
  predict(
    best_model,
    XData = gradient$XDataNew,
    expected = TRUE,
    studyDesign = gradient$studyDesignNew,
    ranLevels = gradient$rLNew)

n_runs <- nChains*samples

pred_df_ff <-
  do.call('rbind', predY_region) %>%
  as_tibble() %>%
  mutate(
    reg = rep(
      gradient$XDataNew$region,
      n_runs),
    run = rep(1:n_runs, each = 4)) %>%
  pivot_longer(
    cols = !c(reg,run),
    values_to = 'occurrence',
    names_to = 'species') %>%
  filter(species %in% annual_10) %>%
  #filter(species %in% perennial_10) %>%
  mutate(species = species %>%
           str_replace('_', ' ')) %>%
  mutate(
    Species_f =
      factor(
        species,
        levels = unique(.$species)))

reg_native <-
  pred_df_ff %>%
  ggplot(
    aes(
      x = reg,
      y = occurrence)) +
  geom_jitter(
    alpha = 0.03,
    aes(group = run),
    key_glyph = 'rect') +
  facet_wrap(
    ~Species_f,
    nrow = 4) +
  xlab('Region') +
  ylab('Probability of occurrence') +
  scale_x_discrete(
    labels = c(
      'cd' = 'Colorado Desert',
      'md' = 'Mojave Desert',
      'scm' = 'Southern California Mountains',
      'sd' = 'Sonoran Desert')) +
  guides(color = guide_legend(
    override.aes = list(alpha = 1))) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 8),
    panel.background = element_rect(fill = NA),
    panel.border = element_rect(
      fill = NA,
      linewidth = 0.75),
    legend.justification = c(1,0),
    legend.position = 'right',
    legend.title = element_blank(),
    strip.text = element_text(
      face = 'italic',
      size = 12,
      family = 'Times'),
    strip.background = element_rect(
      colour = 'black',
      size = 1.5,
      linetype = 'solid',
      fill = 'lightblue'))

reg_native

ggsave(
  reg_native,
  filename = 'output/figures/annual_gradient_reg.png',
  #filename = 'output/figures/perennial_gradient_reg.png',
  width = 20,
  height = 10)

# species co-occurrence ---------------------------------------------------

OmegaCor <-
  computeAssociations(best_model)

# supportLevel <- 0.89
#
# toPlot <-
#   ((OmegaCor[[1]]$support > supportLevel) +
#      (OmegaCor[[1]]$support < (1 - supportLevel)) > 0) *
#   OmegaCor[[1]]$mean
#
# toPlot_p <-
#   ((OmegaCor[[1]]$support > supportLevel) +
#      (OmegaCor[[1]]$support < (1 - supportLevel)) > 0) *
#   OmegaCor[[1]]$support

hmdf_mean <-
  OmegaCor[[1]]$mean %>%
  as.data.frame() %>%
  select(all_of(annual_10))

hmdf_mean <- hmdf_mean[annual_10,]

# as.matrix

omega_plot <-
  ggcorrplot::ggcorrplot(
    hmdf_mean,
    type = 'lower',
    hc.order = TRUE,
    title = 'Occurrence')

ggsave(
  omega_plot,
  # filename = 'output/figures/species_associations_annual.jpg',
  filename = 'output/figures/species_associations_perennial.jpg',
  bg = 'white',
  width = 18,
  height = 18)

# extract species ---------------------------------------------------------

sp <-
  "Achyronychia_cooperi"
   #"Antirrhinum_filipes"
   #"Camissoniopsis_pallida"
   #"Caulanthus_lasiophyllus"
   #"Cistanthe_ambigua"
   #"Eremothera_boothii"
   #"Langloisia_setosissima"
   #"Logfia_depressa"
   #"Nemacladus_glanduliferus"
   #"Oligomeris_linifolia",
   #"Agave_deserti"
   #"Atriplex_hymenelytra"
   #"Cylindropuntia_ganderi"
   #"Ericameria_linearifolia"
   #"Galium_stellatum"
   #"Hesperoyucca_whipplei"
   #"Krameria_bicolor"
   #"Larrea_tridentata"
   #"Penstemon_centranthifolius"
   #"Psorothamnus_schottii"

# species that correlate more than 80%

correlation <-
  OmegaCor[[1]]$mean %>%
  as.data.frame() %>%
  select(all_of(sp)) %>%
  rownames_to_column('species') %>%
  as_tibble() %>%
  rename(correlation = 2) %>%
  filter(correlation > 0.80) %>%
  mutate(
    species =
      str_replace(
        species,
        '_',
        ' ')) %>%
  arrange(species) %>%
  print(n = 100)

# species below 80% correlation

 OmegaCor[[1]]$mean %>%
  as.data.frame() %>%
  select(sp) %>%
  rownames_to_column('species') %>%
  as_tibble() %>%
  select(species)%>%
  mutate(
    species =
      str_replace(
        species,
        '_',
        ' ')) %>%
  anti_join(correlation) %>%
  arrange(species) %>%
  print(n = 100)
