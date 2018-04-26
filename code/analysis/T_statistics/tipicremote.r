# Script to calculate t-stats remotely.
# QDR, NEON ITV, copied from the tipic.r script on 24 Aug 2016

iucn <- read.csv('/mnt/research/neon/external_data/final_external_data/IUCN_mammal_ranges.csv')
load('/mnt/research/neon/final_data/allorganismal2016Aug24.r')
load('/mnt/research/neon/final_data/allsiteplot2016Aug24.r')

# 0. Process by plotID (mostly the same as above except plotID is added as a grouping variable.)

library(cati)
library(dplyr)
library(lubridate)
mam_capture <- mutate(mam_capture, 
                      individualandtag = pmin(as.character(individualID), as.character(tagID), na.rm=TRUE),
                      year = year(date))

mam_noID <- mam_capture[mam_capture$individualandtag=='', c('year', 'siteID','plotID','taxonID','family', 'individualandtag','hindfootLength','earLength','tailLength','totalLength','weight','sex','lifeStage')]

mam_grp <- filter(mam_capture, individualandtag != '') %>%
  group_by(year, siteID, plotID, taxonID, family, individualandtag)

mam_byindiv <- summarize_at(mam_grp, vars(hindfootLength, earLength, tailLength, totalLength, weight), median, na.rm=TRUE)
mam_byindiv_class <- do(mam_grp, data.frame(sex = .$sex[1], lifeStage = names(sort(table(.$lifeStage),decreasing=TRUE))[1]))

mam_byindiv <- cbind(as.data.frame(mam_byindiv), as.data.frame(mam_byindiv_class[,c('sex', 'lifeStage')]))                  
mam_byindiv <-rbind(mam_byindiv, mam_noID)


# Merge covariates with mammal data for plotting
mam_capture_sitemerge <- merge(mam_byindiv, neonsitedata[,-c(3:5)], sort=FALSE, all.x=TRUE)
mam_adults <- subset(mam_capture_sitemerge, lifeStage=='adult')

# Split into years
i2014 <- mam_adults$year == 2014
i2015 <- mam_adults$year == 2015

# Convert adult mammal data to cati format
traits_mam14 <- mam_adults[i2014, c('hindfootLength', 'earLength', 'tailLength', 'totalLength', 'weight')]
traits_mam15 <- mam_adults[i2015, c('hindfootLength', 'earLength', 'tailLength', 'totalLength', 'weight')]



# 1. Define IUCN Regional Pools

mamnames <- unique(mam_capture$scientificName) 
mamnames[mamnames %in% iucn$binomial]
mamnames[!mamnames %in% iucn$binomial] # These include taxa not identified to species as well as two additional species: Tamias alpinus and Rattus rattus. Tamias alpinus should only be in California sites (although it is recorded as a single individual from Utah), and Rattus rattus should be in all sites. Put all the genus-level traits into all the regional species pools.

# Fortunately, only 470 rows in the mam_capture dataframe are not accounted for in IUCN.
#sum(mam_capture$scientificName %in% mamnames[!mamnames %in% iucn$binomial])

nomatchspp <- mamnames[!mamnames %in% iucn$binomial]
nomatchspp <- nomatchspp[c(1:5,7:11,13)] # Only the ones from all sites


sites14 <- unique(mam_adults$siteID[i2014])
sites15 <- unique(mam_adults$siteID[i2015])
plots14 <- unique(mam_adults$plotID[i2014])
plots15 <- unique(mam_adults$plotID[i2015])

reg_pool <- list()

for (plot in unique(mam_capture$plotID)) {
  site <- substr(plot, 1, 4)
  iucn_spp <- as.character(iucn$binomial[iucn[,site]])
  neon_spp <- mam_capture$scientificName[mam_capture$siteID == site]
  regional_spp <- unique(c(iucn_spp, neon_spp, nomatchspp))
  regional_spp_id <- unique(mam_capture$taxonID[mam_capture$scientificName %in% regional_spp])
  reg_pool <- within(reg_pool, assign(plot, value=mam_adults$taxonID %in% regional_spp_id))
}

plotregpoollist_mam14_iucn <- lapply(reg_pool[plots14], function(x) mam_adults[x, c('hindfootLength', 'earLength', 'tailLength', 'totalLength', 'weight')])
plotregpoollist_mam15_iucn <- lapply(reg_pool[plots15], function(x) mam_adults[x, c('hindfootLength', 'earLength', 'tailLength', 'totalLength', 'weight')])

# 2. Define site-only Regional Pools

plotregpoollist_mam14_siteonly <- lapply(plots14, function(x) mam_adults[mam_adults$siteID == substr(x, 1, 4), c('hindfootLength', 'earLength', 'tailLength', 'totalLength', 'weight')])
plotregpoollist_mam15_siteonly <- lapply(plots15, function(x) mam_adults[mam_adults$siteID == substr(x, 1, 4), c('hindfootLength', 'earLength', 'tailLength', 'totalLength', 'weight')])

names(plotregpoollist_mam14_siteonly) <- plots14
names(plotregpoollist_mam15_siteonly) <- plots15

# 3. Calculate T-stats for each year and each regional pool type

tstats_mam_plots14_iucnpool <- Tstats(traits = traits_mam14, ind.plot = factor(mam_adults$plotID[i2014]), sp = factor(mam_adults$taxonID[i2014]), reg.pool = plotregpoollist_mam14_iucn, nperm = 999)
tstats_mam_plots15_iucnpool <- Tstats(traits = traits_mam15, ind.plot = factor(mam_adults$plotID[i2015]), sp = factor(mam_adults$taxonID[i2015]), reg.pool = plotregpoollist_mam15_iucn, nperm = 999)
tstats_mam_plots14_sitepool <- Tstats(traits = traits_mam14, ind.plot = factor(mam_adults$plotID[i2014]), sp = factor(mam_adults$taxonID[i2014]), reg.pool = plotregpoollist_mam14_siteonly, nperm = 999)
tstats_mam_plots15_sitepool <- Tstats(traits = traits_mam15, ind.plot = factor(mam_adults$plotID[i2015]), sp = factor(mam_adults$taxonID[i2015]), reg.pool = plotregpoollist_mam15_siteonly, nperm = 999)

save(tstats_mam_plots14_iucnpool, tstats_mam_plots15_iucnpool, tstats_mam_plots14_sitepool, tstats_mam_plots15_sitepool, file = '~/data/tstatsmambyplot.r')