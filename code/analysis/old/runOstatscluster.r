# Overlap Null Models (without T-stats) run on the mammal data.
# Author: QDR
# Project: NEON ITV
# Created: 30 Sep 2016
# Last modified: 03 Oct 2016

# Set number of null model iterations here:
N_PERM <- 999

# Run on cluster: Load data 

# Load neon data
iucn <- read.csv('/mnt/research/neon/external_data/final_external_data/IUCN_mammal_ranges.csv')
load('/mnt/research/neon/final_data/allorganismal2016Aug24.r')
load('/mnt/research/neon/final_data/allsiteplot2016Aug24.r')

mammalTax <- read.csv('/mnt/research/neon/final_data/mammals/NEON_mam_taxonomy.csv')
mammalTraits <- read.csv('/mnt/research/neon/external_data/final_external_data/NEON_miscmammaltraits.csv')

# Source O-stats calc script
source('~/code/Ostats.r')
# Source community density overlap script
source('~/code/densityoverlap.r')

# Clean mammal data

library(dplyr)
library(lubridate)
mam_capture <- mutate(mam_capture, 
                      individualandtag = pmin(as.character(individualID), as.character(tagID), na.rm=TRUE),
                      year = year(date)) %>%
  left_join(mammalTax[,c('taxonID','order')]) %>% # Added 21 Sep. Gets rid of everything but rodents.
  filter(order == 'Rodentia')

mam_noID <- mam_capture[mam_capture$individualandtag=='', c('year', 'siteID','plotID', 'taxonID','family', 'individualandtag','hindfootLength','earLength','tailLength','totalLength','weight','sex','lifeStage')]

mam_grp <- filter(mam_capture, individualandtag != '') %>%
  group_by(year, siteID, plotID, taxonID, family, individualandtag)

mam_byindiv <- summarize_at(mam_grp, vars(hindfootLength, earLength, tailLength, totalLength, weight), median, na.rm=TRUE)
mam_byindiv_class <- do(mam_grp, data.frame(sex = .$sex[1], lifeStage = names(sort(table(.$lifeStage),decreasing=TRUE))[1]))

mam_byindiv <- cbind(as.data.frame(mam_byindiv), as.data.frame(mam_byindiv_class[,c('sex', 'lifeStage')]))                  
mam_byindiv <-rbind(mam_byindiv, mam_noID)

mam_capture_sitemerge <- merge(mam_byindiv, neonsitedata[,-c(3:5)], sort=FALSE, all.x=TRUE) %>% mutate(logweight = log(weight))

# Split into years
i2014 <- mam_capture_sitemerge$year == 2014
i2015 <- mam_capture_sitemerge$year == 2015
iall <- rep(TRUE, length(i2014))

# Convert adult mammal data to trait-only format
traits_mam14 <- mam_capture_sitemerge[i2014, c('hindfootLength', 'weight', 'logweight')]
traits_mam15 <- mam_capture_sitemerge[i2015, c('hindfootLength', 'weight', 'logweight')]
traits_mamallyears <- mam_capture_sitemerge[, c('hindfootLength', 'weight', 'logweight')]

# Define IUCN Regional Pools

mamnames <- unique(mam_capture$scientificName) 
#mamnames[mamnames %in% iucn$binomial]
#mamnames[!mamnames %in% iucn$binomial] 
# These include taxa not identified to species as well as two additional species: Tamias alpinus and Rattus rattus. Tamias alpinus should only be in California sites (although it is recorded as a single individual from Utah), and Rattus rattus should be in all sites. Put all the genus-level traits into all the regional species pools.

# Fortunately, only 343 (as of Sep 21) rows in the mam_capture dataframe are not accounted for in IUCN.
#sum(mam_capture$scientificName %in% mamnames[!mamnames %in% iucn$binomial])

nomatchspp <- mamnames[!mamnames %in% iucn$binomial]
# Get rid of the individual Tamias alpinus which should not be in regional species pools
nomatchspp <- nomatchspp[!nomatchspp %in% c('Tamias alpinus')]

sites14 <- unique(mam_capture_sitemerge$siteID[i2014])
sites15 <- unique(mam_capture_sitemerge$siteID[i2015])
plots14 <- unique(mam_capture_sitemerge$plotID[i2014])
plots15 <- unique(mam_capture_sitemerge$plotID[i2015])
sitesallyears <- unique(mam_capture_sitemerge$siteID)
plotsallyears <- unique(mam_capture_sitemerge$plotID)

# Get regional pool for each site

reg_pool <- list()

for (site in unique(mam_capture$siteID)) {
  iucn_spp <- as.character(iucn$binomial[iucn[,site]])
  neon_spp <- mam_capture$scientificName[mam_capture$siteID == site]
  regional_spp <- unique(c(iucn_spp, neon_spp, nomatchspp))
  regional_spp_id <- unique(mam_capture$taxonID[mam_capture$scientificName %in% regional_spp])
  reg_pool <- within(reg_pool, assign(site, value=mam_capture_sitemerge$taxonID %in% regional_spp_id))
}

siteregpoollist_mam14_iucn <- lapply(reg_pool[sites14], function(x) mam_capture_sitemerge[x, c('hindfootLength', 'weight', 'logweight')])
siteregpoollist_mam15_iucn <- lapply(reg_pool[sites15], function(x) mam_capture_sitemerge[x, c('hindfootLength', 'weight', 'logweight')])
siteregpoollist_mamallyears_iucn <- lapply(reg_pool, function(x) mam_capture_sitemerge[x, c('hindfootLength', 'weight', 'logweight')])

siteregpoolsp_mam14_iucn <- lapply(reg_pool[sites14], function(x) mam_capture_sitemerge[x, c('taxonID')])
siteregpoolsp_mam15_iucn <- lapply(reg_pool[sites15], function(x) mam_capture_sitemerge[x, c('taxonID')])
siteregpoolsp_mamallyears_iucn <- lapply(reg_pool, function(x) mam_capture_sitemerge[x, c('taxonID')])

# Get task ID
task <- as.numeric(Sys.getenv('PBS_ARRAYID'))

if (task == 1) {
  Ostats_bysite2014 <- Ostats(traits = traits_mam14, plots = factor(mam_capture_sitemerge$siteID[i2014]), sp = factor(mam_capture_sitemerge$taxonID[i2014]), reg_pool_traits = siteregpoollist_mam14_iucn, reg_pool_sp = siteregpoolsp_mam14_iucn, nperm = N_PERM)
  save(Ostats_bysite2014, file = '~/data/Ostats_bysite2014.r')
}
if (task == 2) {
  Ostats_bysite2015 <- Ostats(traits = traits_mam15, plots = factor(mam_capture_sitemerge$siteID[i2015]), sp = factor(mam_capture_sitemerge$taxonID[i2015]), reg_pool_traits = siteregpoollist_mam15_iucn, reg_pool_sp = siteregpoolsp_mam15_iucn, nperm = N_PERM)
  save(Ostats_bysite2015, file = '~/data/Ostats_bysite2015.r')
}
if (task == 3) {
  Ostats_bysiteallyears <- Ostats(traits = traits_mamallyears, plots = factor(mam_capture_sitemerge$siteID[iall]), sp = factor(mam_capture_sitemerge$taxonID[iall]), reg_pool_traits = siteregpoollist_mamallyears_iucn, reg_pool_sp = siteregpoolsp_mamallyears_iucn, nperm = N_PERM)
  save(Ostats_bysiteallyears, file = '~/data/Ostats_bysiteallyears.r')
}
if (task == 4) {
  Ostats_regional2014 <- Ostats_regional(traits = traits_mam14, plots = factor(mam_capture_sitemerge$siteID[i2014]), sp = factor(mam_capture_sitemerge$taxonID[i2014]), reg_pool_traits = siteregpoollist_mam14_iucn, reg_pool_sp = siteregpoolsp_mam14_iucn, nperm = N_PERM)
  save(Ostats_regional2014, file = '~/data/Ostats_regional2014.r')
}
if (task == 5) {
  Ostats_regional2015 <-  Ostats_regional(traits = traits_mam15, plots = factor(mam_capture_sitemerge$siteID[i2015]), sp = factor(mam_capture_sitemerge$taxonID[i2015]), reg_pool_traits = siteregpoollist_mam15_iucn, reg_pool_sp = siteregpoolsp_mam15_iucn, nperm = N_PERM)
  save(Ostats_regional2015, file = '~/data/Ostats_regional2015.r')
}
if (task == 6) {
  Ostats_regionalallyears <- Ostats_regional(traits = traits_mamallyears, plots = factor(mam_capture_sitemerge$siteID[iall]), sp = factor(mam_capture_sitemerge$taxonID[iall]), reg_pool_traits = siteregpoollist_mamallyears_iucn, reg_pool_sp = siteregpoolsp_mamallyears_iucn, nperm = N_PERM)
  save(Ostats_regionalallyears, file = '~/data/Ostats_regionalallyears.r')
}
