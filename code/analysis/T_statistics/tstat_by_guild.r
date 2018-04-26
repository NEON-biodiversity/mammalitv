# T-stats by guild
# Author: QDR
# Project: NEON ITV
# Created: 21 Sep 2016
# Last modified: 22 Sep 2016

# Modified 22 Sep: The t-stats did not calculate for the previous code (site-only RPs). I got rid of all the plots that don't have enough individuals to get a meaningful t-statistic.
# Also add in the individuals' guilds that are only identified to species.

# Modified version of the tipic.r script that is split by guild, at the subplot level.
# Also set up to run remotely.

library(cati)

### LOAD DATA ON CLUSTER

# Load neon data
iucn <- read.csv('/mnt/research/neon/external_data/final_external_data/IUCN_mammal_ranges.csv')
load('/mnt/research/neon/final_data/allorganismal2016Aug24.r')
load('/mnt/research/neon/final_data/allsiteplot2016Aug24.r')

mammalTax <- read.csv('/mnt/research/neon/final_data/mammals/NEON_mam_taxonomy.csv')
mammalTraits <- read.csv('/mnt/research/neon/external_data/final_external_data/NEON_miscmammaltraits.csv')

### LOAD DATA LOCALLY

# iucn <- read.csv('C:/Users/Q/Dropbox/neon/data/external_datasets/IUCN_mammal_ranges.csv')
# load('allorganismal_latest.r')
# load('allsiteplot_latest.r')
# 
# mammalTax <- read.csv('C:/Users/Q/Dropbox/neon/protocols/taxonomy/NEON_mam_taxonomy.csv')
# mammalTraits <- read.csv('C:/Users/Q/Dropbox/neon/data/external_datasets/NEON_miscmammaltraits.csv')

# Clean mammal data to get just adults.
# Take the median value of the traits for each individual if it was recaptured
library(dplyr)
library(lubridate)
mam_capture <- mutate(mam_capture, 
                      individualandtag = pmin(as.character(individualID), as.character(tagID), na.rm=TRUE),
                      year = year(date)) %>%
  left_join(mammalTax[,c('taxonID','order')]) %>% # Added 21 Sep. Gets rid of everything but rodents.
  filter(order == 'Rodentia')

mammalGuilds <- mammalTraits %>% 
  mutate(scientificName = paste(Genus, Species)) %>%
  select(scientificName, Pineda_Main_food)

mammalGuilds <- rbind(mammalGuilds, data.frame(scientificName=c("Dipodomys sp.", "Glaucomys sp.", "Microtus sp.", "Neotoma sp.", 
                                                                "Perognathus sp.", "Peromyscus sp.", "Reithrodontomys sp."),
                                               Pineda_Main_food=c('Generalist','Generalist','Herbivore','Herbivore','Generalist','Granivore','Generalist')))

mam_capture <- left_join(mam_capture, mammalGuilds[-c(89,118,62,63),], by = 'scientificName')

mam_noID <- mam_capture[mam_capture$individualandtag=='', c('year', 'siteID','plotID', 'taxonID','family', 'individualandtag','hindfootLength','earLength','tailLength','totalLength','weight','sex','lifeStage','Pineda_Main_food')]

mam_grp <- filter(mam_capture, individualandtag != '') %>%
  group_by(year, siteID, plotID, taxonID, family, individualandtag)

mam_byindiv <- summarize_at(mam_grp, vars(hindfootLength, earLength, tailLength, totalLength, weight), median, na.rm=TRUE)
mam_byindiv_class <- do(mam_grp, data.frame(sex = .$sex[1], lifeStage = names(sort(table(.$lifeStage),decreasing=TRUE))[1], Pineda_Main_food=.$Pineda_Main_food[1]))

mam_byindiv <- cbind(as.data.frame(mam_byindiv), as.data.frame(mam_byindiv_class[,c('sex', 'lifeStage','Pineda_Main_food')]))                  
mam_byindiv <-rbind(mam_byindiv, mam_noID)

mam_capture_sitemerge <- merge(mam_byindiv, neonsitedata[,-c(3:5)], sort=FALSE, all.x=TRUE) %>% mutate(logweight = log(weight))

# Split into years
i2014 <- mam_capture_sitemerge$year == 2014
i2015 <- mam_capture_sitemerge$year == 2015

# Split into functional groupings for 2015 only
i2015generalistsgranivores <- with(mam_capture_sitemerge, year == 2015 & Pineda_Main_food %in% c('Generalist','Granivore'))
i2015generalistsherbivores <- with(mam_capture_sitemerge, year == 2015 & Pineda_Main_food %in% c('Generalist','Herbivore'))
i2015threeguilds <- with(mam_capture_sitemerge, year == 2015 & Pineda_Main_food %in% c('Generalist','Granivore','Herbivore'))

# (amended Sep22): Get rid of all plots with 4 or fewer individuals.

plottable_gg <- table(factor(mam_capture_sitemerge$plotID[i2015generalistsgranivores]))
plottable_gh <- table(factor(mam_capture_sitemerge$plotID[i2015generalistsherbivores]))
plottable_ggh <- table(factor(mam_capture_sitemerge$plotID[i2015threeguilds]))

plotuse_gg <- names(plottable_gg)[plottable_gg > 4]
plotuse_gh <- names(plottable_gh)[plottable_gh > 4]
plotuse_ggh <- names(plottable_ggh)[plottable_ggh > 4]

i2015generalistsgranivores <- i2015generalistsgranivores & mam_capture_sitemerge$plotID %in% plotuse_gg
i2015generalistsherbivores <- i2015generalistsherbivores & mam_capture_sitemerge$plotID %in% plotuse_gh
i2015threeguilds <- i2015threeguilds & mam_capture_sitemerge$plotID %in% plotuse_ggh

# Convert adult mammal data to cati format
traits_mam14 <- mam_capture_sitemerge[i2014, c('hindfootLength', 'earLength', 'tailLength', 'totalLength', 'weight', 'logweight')]
traits_mam15 <- mam_capture_sitemerge[i2015, c('hindfootLength', 'earLength', 'tailLength', 'totalLength', 'weight', 'logweight')]
traits_mam15gg <- mam_capture_sitemerge[i2015generalistsgranivores, c('hindfootLength', 'earLength', 'tailLength', 'totalLength', 'weight', 'logweight')]
traits_mam15gh <- mam_capture_sitemerge[i2015generalistsherbivores, c('hindfootLength', 'earLength', 'tailLength', 'totalLength', 'weight', 'logweight')]
traits_mam15ggh <- mam_capture_sitemerge[i2015threeguilds, c('hindfootLength', 'earLength', 'tailLength', 'totalLength', 'weight', 'logweight')]


# 1. Define IUCN Regional Pools

mamnames <- unique(mam_capture$scientificName) 
mamnames[mamnames %in% iucn$binomial]
mamnames[!mamnames %in% iucn$binomial] # These include taxa not identified to species as well as two additional species: Tamias alpinus and Rattus rattus. Tamias alpinus should only be in California sites (although it is recorded as a single individual from Utah), and Rattus rattus should be in all sites. Put all the genus-level traits into all the regional species pools.

# Fortunately, only 343 (as of Sep 21) rows in the mam_capture dataframe are not accounted for in IUCN.
sum(mam_capture$scientificName %in% mamnames[!mamnames %in% iucn$binomial])

nomatchspp <- mamnames[!mamnames %in% iucn$binomial]
# Get rid of the individual Tamias alpinus which should not be in regional species pools
nomatchspp <- nomatchspp[!nomatchspp %in% c('Tamias alpinus')]

sites14 <- unique(mam_capture_sitemerge$siteID[i2014])
sites15 <- unique(mam_capture_sitemerge$siteID[i2015])
plots14 <- unique(mam_capture_sitemerge$plotID[i2014])
plots15 <- unique(mam_capture_sitemerge$plotID[i2015])

sites15gg <- unique(mam_capture_sitemerge$siteID[i2015generalistsgranivores])
plots15gg <- unique(mam_capture_sitemerge$plotID[i2015generalistsgranivores])
sites15gh <- unique(mam_capture_sitemerge$siteID[i2015generalistsherbivores])
plots15gh <- unique(mam_capture_sitemerge$plotID[i2015generalistsherbivores])
sites15ggh <- unique(mam_capture_sitemerge$siteID[i2015threeguilds])
plots15ggh <- unique(mam_capture_sitemerge$plotID[i2015threeguilds])

reg_pool <- list()

for (plot in unique(mam_capture$plotID)) {
  site <- substr(plot, 1, 4)
  iucn_spp <- as.character(iucn$binomial[iucn[,site]])
  neon_spp <- mam_capture$scientificName[mam_capture$siteID == site]
  regional_spp <- unique(c(iucn_spp, neon_spp, nomatchspp))
  regional_spp_id <- unique(mam_capture$taxonID[mam_capture$scientificName %in% regional_spp])
  reg_pool <- within(reg_pool, assign(plot, value=mam_capture_sitemerge$taxonID %in% regional_spp_id))
}

plotregpoollist_mam14_iucn <- lapply(reg_pool[plots14], function(x) mam_capture_sitemerge[x, c('hindfootLength', 'earLength', 'tailLength', 'totalLength', 'weight', 'logweight')])
plotregpoollist_mam15_iucn <- lapply(reg_pool[plots15], function(x) mam_capture_sitemerge[x, c('hindfootLength', 'earLength', 'tailLength', 'totalLength', 'weight', 'logweight')])
plotregpoollist_mam15gg_iucn <- lapply(reg_pool[plots15gg], function(x) mam_capture_sitemerge[x, c('hindfootLength', 'earLength', 'tailLength', 'totalLength', 'weight', 'logweight')])
plotregpoollist_mam15gh_iucn <- lapply(reg_pool[plots15gh], function(x) mam_capture_sitemerge[x, c('hindfootLength', 'earLength', 'tailLength', 'totalLength', 'weight', 'logweight')])
plotregpoollist_mam15ggh_iucn <- lapply(reg_pool[plots15ggh], function(x) mam_capture_sitemerge[x, c('hindfootLength', 'earLength', 'tailLength', 'totalLength', 'weight', 'logweight')])


# 2. Define site-only Regional Pools

plotregpoollist_mam14_siteonly <- lapply(plots14, function(x) {
  spp <- unique(mam_capture_sitemerge$taxonID[mam_capture_sitemerge$siteID == substr(x, 1, 4)])
  mam_capture_sitemerge[mam_capture_sitemerge$taxonID %in% spp, c('hindfootLength', 'earLength', 'tailLength', 'totalLength', 'weight', 'logweight')]
})
plotregpoollist_mam15_siteonly <- lapply(plots15, function(x) {
  spp <- unique(mam_capture_sitemerge$taxonID[mam_capture_sitemerge$siteID == substr(x, 1, 4)])
  mam_capture_sitemerge[mam_capture_sitemerge$taxonID %in% spp, c('hindfootLength', 'earLength', 'tailLength', 'totalLength', 'weight', 'logweight')]
})
plotregpoollist_mam15gg_siteonly <- lapply(plots15gg, function(x) {
  spp <- unique(mam_capture_sitemerge$taxonID[mam_capture_sitemerge$siteID == substr(x, 1, 4)])
  mam_capture_sitemerge[mam_capture_sitemerge$taxonID %in% spp, c('hindfootLength', 'earLength', 'tailLength', 'totalLength', 'weight', 'logweight')]
})
plotregpoollist_mam15gh_siteonly <- lapply(plots15gh, function(x) {
  spp <- unique(mam_capture_sitemerge$taxonID[mam_capture_sitemerge$siteID == substr(x, 1, 4)])
  mam_capture_sitemerge[mam_capture_sitemerge$taxonID %in% spp, c('hindfootLength', 'earLength', 'tailLength', 'totalLength', 'weight', 'logweight')]
})
plotregpoollist_mam15ggh_siteonly <- lapply(plots15ggh, function(x) {
  spp <- unique(mam_capture_sitemerge$taxonID[mam_capture_sitemerge$siteID == substr(x, 1, 4)])
  mam_capture_sitemerge[mam_capture_sitemerge$taxonID %in% spp, c('hindfootLength', 'earLength', 'tailLength', 'totalLength', 'weight', 'logweight')]
})

names(plotregpoollist_mam14_siteonly) <- plots14
names(plotregpoollist_mam15_siteonly) <- plots15

names(plotregpoollist_mam15gg_siteonly) <- plots15gg
names(plotregpoollist_mam15gh_siteonly) <- plots15gh
names(plotregpoollist_mam15ggh_siteonly) <- plots15ggh


# 3. Calculate T-stats for  each regional pool type and each guild grouping (2015 only)

task <- as.numeric(Sys.getenv('PBS_ARRAYID'))
print(task)

traits_i <- list(traits_mam15gg, traits_mam15gh, traits_mam15ggh, traits_mam15gg, traits_mam15gh, traits_mam15ggh)[[task]]
i_i <- list(i2015generalistsgranivores, i2015generalistsherbivores, i2015threeguilds, i2015generalistsgranivores, i2015generalistsherbivores, i2015threeguilds)[[task]]
pool_i <- list(plotregpoollist_mam15gg_iucn, plotregpoollist_mam15gh_iucn, plotregpoollist_mam15ggh_iucn, plotregpoollist_mam15gg_siteonly, plotregpoollist_mam15gh_siteonly, plotregpoollist_mam15ggh_siteonly)[[task]]
name_i <- c('tstats_2015_gg_iucn', 'tstats_2015_gh_iucn', 'tstats_2015_ggh_iucn', 'tstats_2015_gg_siteonly', 'tstats_2015_gh_siteonly','tstats_2015_ggh_siteonly')[task]

tstats_i <- Tstats(traits = traits_i, ind.plot = factor(mam_capture_sitemerge$plotID[i_i]), sp = factor(mam_capture_sitemerge$taxonID[i_i]), reg.pool = pool_i, nperm = 999)
assign(name_i, value = tstats_i)

save(list = name_i, file = paste0('~/data/',name_i,'.r'))