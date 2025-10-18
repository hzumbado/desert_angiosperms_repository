# project structure

folders <-
  c('data/raw',
    'data/processed',
    'shapefiles',
    'rasters',
    'scripts/richness',
    'scripts/jsdms',
    'scripts/sdms',
    'output/figures',
    'output/jsdms/',
    'output/sdms')

sapply(
  folders,
  FUN = dir.create,
  recursive = TRUE)
