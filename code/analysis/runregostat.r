source('code/setupostats.r')

 Ostats_regional2014 <- Ostats_regional(traits = traits_mam14, plots = factor(mam_capture_sitemerge$siteID[i2014]), sp = factor(mam_capture_sitemerge$taxonID[i2014]), reg_pool_traits = siteregpoollist_mam14_iucn, reg_pool_sp = siteregpoolsp_mam14_iucn, nperm = N_PERM)
 Ostats_regional2015 <-  Ostats_regional(traits = traits_mam15, plots = factor(mam_capture_sitemerge$siteID[i2015]), sp = factor(mam_capture_sitemerge$taxonID[i2015]), reg_pool_traits = siteregpoollist_mam15_iucn, reg_pool_sp = siteregpoolsp_mam15_iucn, nperm = N_PERM)
 Ostats_regionalallyears <- Ostats_regional(traits = traits_mamallyears, plots = factor(mam_capture_sitemerge$siteID[iall]), sp = factor(mam_capture_sitemerge$taxonID[iall]), reg_pool_traits = siteregpoollist_mamallyears_iucn, reg_pool_sp = siteregpoolsp_mamallyears_iucn, nperm = N_PERM)

 Ostats_regcontinent2014 <- Ostats_regional(traits = traits_mam14, plots = factor(mam_capture_sitemerge$siteID[i2014]), sp = factor(mam_capture_sitemerge$taxonID[i2014]), reg_pool_traits = continentpool14, reg_pool_sp = continentsp14, nperm = N_PERM)
 Ostats_regcontinent2015 <-  Ostats_regional(traits = traits_mam15, plots = factor(mam_capture_sitemerge$siteID[i2015]), sp = factor(mam_capture_sitemerge$taxonID[i2015]), reg_pool_traits = continentpool15, reg_pool_sp = continentsp15, nperm = N_PERM)
 Ostats_regcontinentallyears <- Ostats_regional(traits = traits_mamallyears, plots = factor(mam_capture_sitemerge$siteID[iall]), sp = factor(mam_capture_sitemerge$taxonID[iall]), reg_pool_traits = continentpoolallyears, reg_pool_sp = continentspallyears, nperm = N_PERM)

save(Ostats_regcontinent2014, Ostats_regcontinent2015, Ostats_regcontinentallyears, Ostats_regional2014, Ostats_regional2015, Ostats_regionalallyears, file = '~/data/Ostats_allregional.r')

