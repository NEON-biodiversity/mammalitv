# Extract richness for NEON sites from NEON and from the Ellis map.
# For ecography revision.
# QDR 08 Nov 2017 (NEON ITV)

# Edit 25 Apr 2018: Change path to HPCC, not Dropbox
data_path <- '/mnt/research/neon'

library(dplyr)

# NEON --------------------------------------------------------------------

# Use slick new API to get all the plant data from NEON.
# Run on 08 Nov 2017.
source('~/GitHub/NEON/code/data_extraction/datapull_neonapi_fns.r')
plantppc_code <- 'DP1.10058.001'
display_neon_filenames(plantppc_code)
plantppc <- pull_all_neon_data(plantppc_code, '100m2Data') # 100 m2 plots.
# Note that there is no longer a different file for the 400 m2 plots. You just have to sum up the presences in the 100 m2 plots.

r100m2 <- plantppc %>%
  group_by(siteID) %>%
  summarize(richness_NEON = length(unique(scientificName)))

# Ellis -------------------------------------------------------------------

library(maptools)
ellis <- readShapeSpatial(fn = file.path(data_path, 'external_data/raw_external_data/other/ellis_2012_shapefile/ellis_2012_l8_dataset_2012_01_17'))

# Get richness for the neon sites, along with their coordinates, and also compare this to the background richness.
load(file.path(data_path, 'final_data/allorganismal_latest.r'))
load(file.path(data_path, 'final_data/allsiteplot_latest.r'))

# Background richness from ellis
neon_ellis <- extract(x = ellis, y = with(neonsitedata, data.frame(x = decimalLongitude, y = decimalLatitude)))


# Join and export ---------------------------------------------------------


plantrichness <- neonsitedata %>% 
  cbind(neon_ellis) %>%
  left_join(r100m2) %>%
  dplyr::select(siteID, richness_NEON, N) %>%
  rename(richness_Ellis = N)

write.csv(plantrichness, file = file.path(data_path, 'external_data/final_external_data/NEON_background_plant_richness.csv'), row.names = FALSE)
