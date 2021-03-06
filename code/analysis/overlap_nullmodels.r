# T-stats and overlap stats run on the cluster
# Author: QDR
# Project: NEON ITV
# Created: 21 Sep 2016
# Last modified: 29 Sep 2016

# Modified 29 Sep (saved as a new script): Added the overlap stat to the null model. Run for each year (2014 and 2015), and all years together. Don't split by guild.
# Modified 26 Sep (saved as a new script): Added my tweaked t-stats function that uses weighted means. Also added each guild separately. Finally, just do the IUCN regional pools, not caring about siteonly for now.
# Modified 22 Sep: The t-stats did not calculate for the previous code (site-only RPs). I got rid of all the plots that don't have enough individuals to get a meaningful t-statistic.
# Also add in the individuals' guilds that are only identified to species.

# Modified version of the tipic.r script that is split by guild, at the subplot level.
# Also set up to run remotely.

library(cati)

N_PERM <- 999 # From now on, edit permutation number here.

### LOAD DATA ON CLUSTER

# Load neon data
iucn <- read.csv('/mnt/research/neon/external_data/final_external_data/IUCN_mammal_ranges.csv')
load('/mnt/research/neon/final_data/allorganismal2016Aug24.r')
load('/mnt/research/neon/final_data/allsiteplot2016Aug24.r')

mammalTax <- read.csv('/mnt/research/neon/final_data/mammals/NEON_mam_taxonomy.csv')
mammalTraits <- read.csv('/mnt/research/neon/external_data/final_external_data/NEON_miscmammaltraits.csv')

# Source tweaked t-stats calc script
source('code/analysis/T_statistics/Tstatsedited.r')
# Source community density overlap script
source('code/analysis/densityoverlap.r')

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
  dplyr::select(scientificName, Pineda_Main_food)

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
iall <- rep(TRUE, length(i2014))

# Convert adult mammal data to cati format
traits_mam14 <- mam_capture_sitemerge[i2014, c('hindfootLength', 'earLength', 'tailLength', 'totalLength', 'weight', 'logweight')]
traits_mam15 <- mam_capture_sitemerge[i2015, c('hindfootLength', 'earLength', 'tailLength', 'totalLength', 'weight', 'logweight')]
traits_mamallyears <- mam_capture_sitemerge[, c('hindfootLength', 'earLength', 'tailLength', 'totalLength', 'weight', 'logweight')]


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
sitesallyears <- unique(mam_capture_sitemerge$siteID)
plotsallyears <- unique(mam_capture_sitemerge$plotID)

# Regional pool definition for each plot

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
plotregpoollist_mamallyears_iucn <- lapply(reg_pool, function(x) mam_capture_sitemerge[x, c('hindfootLength', 'earLength', 'tailLength', 'totalLength', 'weight', 'logweight')])

# Regional pool definition for each site

reg_pool <- list()

for (site in unique(mam_capture$siteID)) {
  iucn_spp <- as.character(iucn$binomial[iucn[,site]])
  neon_spp <- mam_capture$scientificName[mam_capture$siteID == site]
  regional_spp <- unique(c(iucn_spp, neon_spp, nomatchspp))
  regional_spp_id <- unique(mam_capture$taxonID[mam_capture$scientificName %in% regional_spp])
  reg_pool <- within(reg_pool, assign(site, value=mam_capture_sitemerge$taxonID %in% regional_spp_id))
}

siteregpoollist_mam14_iucn <- lapply(reg_pool[sites14], function(x) mam_capture_sitemerge[x, c('hindfootLength', 'earLength', 'tailLength', 'totalLength', 'weight', 'logweight')])
siteregpoollist_mam15_iucn <- lapply(reg_pool[sites15], function(x) mam_capture_sitemerge[x, c('hindfootLength', 'earLength', 'tailLength', 'totalLength', 'weight', 'logweight')])
siteregpoollist_mamallyears_iucn <- lapply(reg_pool, function(x) mam_capture_sitemerge[x, c('hindfootLength', 'earLength', 'tailLength', 'totalLength', 'weight', 'logweight')])



# 3. Calculate T-stats and overlap stats for each year
task <- as.numeric(Sys.getenv('PBS_ARRAYID'))
print(task) # from 1 to 6

traits_i <- list(traits_mam14, traits_mam15, traits_mamallyears, traits_mam14, traits_mam15, traits_mamallyears)[[task]]
indplot_i <- list(factor(mam_capture_sitemerge$siteID[i2014]), factor(mam_capture_sitemerge$siteID[i2015]), factor(mam_capture_sitemerge$siteID[iall]), factor(mam_capture_sitemerge$plotID[i2014]), factor(mam_capture_sitemerge$plotID[i2015]), factor(mam_capture_sitemerge$plotID[iall]))[[task]]
i_i <- list(i2014, i2015, iall, i2014, i2015, iall)[[task]]
pool_i <- list(siteregpoollist_mam14_iucn, siteregpoollist_mam15_iucn, siteregpoollist_mamallyears_iucn, plotregpoollist_mam14_iucn, plotregpoollist_mam15_iucn, plotregpoollist_mamallyears_iucn)[[task]]
name_i <- c('overlap_2014_bysite', 'overlap_2015_bysite', 'overlap_allyears_bysite', 'overlap_2014_byplot', 'overlap_2015_byplot', 'overlap_allyears_byplot')[task]

tstats_i <- QTstats(traits = traits_i, ind.plot = indplot_i, sp = factor(mam_capture_sitemerge$taxonID[i_i]), reg.pool = pool_i, nperm = N_PERM) # It's QTstats because of the weighted mean
assign(name_i, value = tstats_i)

save(list = name_i, file = paste0('~/data/',name_i,'.r'))