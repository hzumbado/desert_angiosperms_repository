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
    'scripts/03-sdm',
    'output/figures',
    'output/models/jsdm',
    'output/models/sdm',
    'output/tables/jsdm',
    'output/tables/sdm')

sapply(
  folders,
  FUN = dir.create,
  recursive = TRUE)
