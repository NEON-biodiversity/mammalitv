# Export data for making figure 3
# After running newloadplotdat.r

# o2015 is already joined with neonsitedata

ovl_dat <- o2015 %>% filter(trait %in% 'logweight') %>% dplyr::select(ostat_norm, siteID)
indiv_dat <- filter(mam_capture_sitemerge, year == 2015) %>% dplyr::select(siteID, taxonID, weight)

library(readr)

write_csv(ovl_dat, '~/Documents/GitHub/NEON_repos/Ostats/internal_docs/overlap_dat.csv')
write_csv(indiv_dat, '~/Documents/GitHub/NEON_repos/Ostats/internal_docs/indiv_dat.csv')
