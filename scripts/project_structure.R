# project structure

folders <-
  c('data/raw',
    'data/processed/richness',
    'data/processed/jsdm',
    'data/processed/sdm',
    'shapefiles',
    'rasters',
    'scripts/01-richness',
    'scripts/02-jsdm',
    'scripts/03-sdm/preprocessing',
    'scripts/03-sdm/modeling',
    'output/figures',
    'output/models/jsdm',
    'output/models/sdm/files/',
    'output/models/sdm/predictions/',
    'output/models/sdm/predictions/average/',
    'output/models/sdm/binary/present/polygons/',
    'output/models/sdm/binary/present/clamped_preds/',
    'output/tables/jsdm',
    'output/tables/sdm')

sapply(
  folders,
  FUN = dir.create,
  recursive = TRUE)
