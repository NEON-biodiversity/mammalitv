
# Initial setup of data for supplement (not to be included in final script)
rm(list=ls())
source('code/vis/newloadplotdat_nullswap.r')


modisdat <- read.csv('C:/Users/Q/Dropbox/neon/data/external_datasets/NEON_modisyearly_withcv.csv', stringsAsFactors = FALSE)
modis_mean <- modisdat %>% group_by(siteID) %>% summarize_at(-(1:3), mean, na.rm=TRUE)

site_covariates <- o2015 %>% 
  filter(trait == 'logweight') %>%
  dplyr::select(1, 25:50, 53, 102:123)

heterodf <- read.csv('C:/Users/Q/Dropbox/neon/data/heterogeneity.csv', stringsAsFactors = FALSE)

site_covariates <- site_covariates %>%
  dplyr::select(-mpd_z, -mntd_z) %>%
  left_join(modis_mean) %>%
  left_join(heterodf)

write.csv(site_covariates, file = 'C:/Users/Q/google_drive/NEON_EAGER/Manuscript/revision/supplemental/site_covariates.csv', row.names = FALSE)
write.csv(mam_capture, file = 'C:/Users/Q/google_drive/NEON_EAGER/Manuscript/revision/supplemental/raw_NEON_mammal_data.csv', row.names = FALSE)
write.csv(mammalTraits, file = 'C:/Users/Q/google_drive/NEON_EAGER/Manuscript/revision/supplemental/mammal_traits.csv', row.names = FALSE)

# Get relevant rows from the mam_capture_sitemerge dataframe
mammal_data_final <- mam_capture_sitemerge %>%
  mutate(logweight = log(weight)) %>%
  filter(year == 2015) %>%
  dplyr::select(domainName, domainID, siteName, siteID, plotID, decimalLatitude, decimalLongitude, taxonID, family, individualandtag, hindfootLength, earLength, tailLength, totalLength, weight, logweight, sex, lifeStage, Pineda_Main_food)

write.csv(mammal_data_final, file = 'C:/Users/Q/google_drive/NEON_EAGER/Manuscript/revision/supplemental/final_NEON_mammal_data.csv', row.names = FALSE)
