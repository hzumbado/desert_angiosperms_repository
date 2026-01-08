# project structure

folders <-
  c('data/raw/richness',
    'data/raw/jsdm',
    'data/raw/sdm',
    'data/processed/richness',
    'data/processed/jsdm',
    'data/processed/sdm',
    'data/processed/sdm/premodel',
    'data/processed/sdm/model',
    'shapefiles',
    'rasters',
    'scripts/01-richness',
    'scripts/02-jsdm',
    'scripts/03-sdm/preprocessing',
    'scripts/03-sdm/modeling',
    'output/figures',
    'output/models/jsdm',
    'output/models/sdm/files',
    'output/models/sdm/predictions/average',
    'output/models/sdm/predictions/present',
    'output/models/sdm/predictions/future',
    'output/models/sdm/binary/present/polygons',
    'output/models/sdm/binary/present/clamped_preds',
    'output/models/sdm/binary/future/polygons',
    'output/models/sdm/binary/future/clamped_preds',
    'output/models/sdm/binary/future/ssp/1_2.6',
    'output/models/sdm/binary/future/ssp/2_4.5',
    'output/models/sdm/binary/future/ssp/5_8.5',
    'output/tables/jsdm',
    'output/tables/sdm')

sapply(
  folders,
  FUN = dir.create,
  recursive = TRUE)
