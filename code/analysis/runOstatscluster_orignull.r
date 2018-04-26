# Overlap Null Models (without T-stats) run on the mammal data. Run on cluster and split by guild.
# Author: QDR
# Project: NEON ITV
# Created: 11 Oct 2016 *from older script

# 13 Jan: Update to include new mam_capture data. Not 2016 tho.
# 09 Jan: Remove opportunistic and non-target species. Add option to just do the weight-shuffling null model.
# 20 Dec: Change weighting statistic to harmonic mean of the abundance pairs, rather than geometric mean.
# 14 Dec: Add three levels of region for the null models. Also don't use anything but 2015 and generalists.
# 25 Oct: Added code to get rid of the insectivorous mammals and the individuals not identified to species.
# 13 Oct: Add weighted median.
# 11 Oct: change to run remotely and split by guild. 
# Modified 10 Oct: Add regional O-stats calculation with continent as regional pool.

# Set number of null model iterations here and which stat is used. Also whether the weights should be shuffled.
STAT <- 'harmonic'
N_PERM <- 999
SW <- FALSE 

# Run on cluster: Load data 

# Load neon data
iucn <- read.csv('/mnt/research/neon/external_data/final_external_data/IUCN_mammal_ranges.csv')
load('/mnt/research/neon/final_data/allorganismal2017Jan13.r')
load('/mnt/research/neon/final_data/allsiteplot2016Aug24.r')

mammalTax <- read.csv('/mnt/research/neon/final_data/mammals/NEON_mam_taxonomy.csv')
mammalTraits <- read.csv('/mnt/research/neon/external_data/final_external_data/NEON_miscmammaltraits.csv')

# Source O-stats calc script
source('~/code/ostatcode/Ostats.r')
# Source community density overlap script
source('~/code/ostatcode/densityoverlap.r')

# Clean mammal data

library(dplyr)
library(lubridate)
mam_capture <- mutate(mam_capture, 
                      individualandtag = pmin(as.character(individualID), as.character(tagID), na.rm=TRUE),
                      year = year(date)) %>%
  left_join(mammalTax[,c('taxonID','order','taxonProtocolCategory')]) %>% # Added 21 Sep. Gets rid of everything but rodents. 09 Jan: Gets rid of opportunistic and non-target species too.
  filter(order == 'Rodentia', taxonProtocolCategory == 'target', year == 2015)

#mammalGuilds <- mammalTraits %>% 
#  mutate(scientificName = paste(Genus, Species)) %>%
#  dplyr::select(scientificName, Pineda_Main_food)
#
#mammalGuilds <- rbind(mammalGuilds, data.frame(scientificName=c("Dipodomys sp.", "Glaucomys sp.", "Microtus sp.", "Neotoma sp.", 
#                                                                "Perognathus sp.", "Peromyscus sp.", "Reithrodontomys sp."),
#                                               Pineda_Main_food=c('Generalist','Generalist','Herbivore','Herbivore','Generalist','Granivore','Generalist')))


#mam_capture <- left_join(mam_capture, mammalGuilds[-c(89,118,62,63),], by = 'scientificName')

mam_noID <- mam_capture[mam_capture$individualandtag=='', c('year', 'siteID','plotID', 'taxonID','family', 'individualandtag','hindfootLength','earLength','tailLength','totalLength','weight','sex','lifeStage')]

mam_grp <- filter(mam_capture, individualandtag != '') %>%
  group_by(year, siteID, plotID, taxonID, family, individualandtag)

mam_byindiv <- summarize_at(mam_grp, vars(hindfootLength, earLength, tailLength, totalLength, weight), median, na.rm=TRUE)
mam_byindiv_class <- do(mam_grp, data.frame(sex = .$sex[1], lifeStage = names(sort(table(.$lifeStage),decreasing=TRUE))[1]))

mam_byindiv <- cbind(as.data.frame(mam_byindiv), as.data.frame(mam_byindiv_class[,c('sex', 'lifeStage')]))                  
mam_byindiv <-rbind(mam_byindiv, mam_noID)

mam_capture_sitemerge <- merge(mam_byindiv, neonsitedata[,-c(3:5)], sort=FALSE, all.x=TRUE) %>% mutate(logweight = log(weight))

# Here, get rid of the insectivores and the mammals not identified to species.
nottospecies <- c('DPSP','GLSP','MISP','NESP','PESP','PGSP','RESP')
mam_capture_sitemerge <- mam_capture_sitemerge %>% filter(!(taxonID %in% nottospecies) &!(siteID %in% c('DSNY','DELA','HEAL'))) # Gets rid of ~700.

# Split into year+guild combination
#i2015 <- mam_capture_sitemerge$year == 2015

# Convert adult mammal data to trait-only format
traits_mam15 <- mam_capture_sitemerge[, c('hindfootLength', 'weight', 'logweight')]

# Define IUCN Regional Pools

mamnames <- unique(mam_capture$scientificName) 

nomatchspp <- mamnames[!mamnames %in% iucn$binomial]
# Get rid of the individual Tamias alpinus which should not be in regional species pools
# nomatchspp <- nomatchspp[!nomatchspp %in% c('Tamias alpinus')]

sites15 <- unique(mam_capture_sitemerge$siteID)
plots15 <- unique(mam_capture_sitemerge$plotID)

# Get local pool for each site

local_pool <- list()

for (site in unique(mam_capture$siteID)) {
  neon_spp <- mam_capture$scientificName[mam_capture$siteID == site]
  regional_spp <- unique(c(neon_spp, nomatchspp))
  regional_spp_id <- unique(mam_capture$taxonID[mam_capture$scientificName %in% regional_spp])
  local_pool <- within(local_pool, assign(site, value=mam_capture_sitemerge$taxonID %in% regional_spp_id))
}

siteregpoollist_mam15_local <- lapply(local_pool[sites15], function(x) mam_capture_sitemerge[x, c('hindfootLength', 'weight', 'logweight')])
siteregpoolsp_mam15_local <- lapply(local_pool[sites15], function(x) mam_capture_sitemerge[x, c('taxonID')])


# Get IUCN regional pool for each site

reg_pool <- list()

for (site in unique(mam_capture$siteID)) {
  iucn_spp <- as.character(iucn$binomial[iucn[,site]])
  neon_spp <- mam_capture$scientificName[mam_capture$siteID == site]
  regional_spp <- unique(c(iucn_spp, neon_spp, nomatchspp))
  regional_spp_id <- unique(mam_capture$taxonID[mam_capture$scientificName %in% regional_spp])
  reg_pool <- within(reg_pool, assign(site, value=mam_capture_sitemerge$taxonID %in% regional_spp_id))
}

siteregpoollist_mam15_iucn <- lapply(reg_pool[sites15], function(x) mam_capture_sitemerge[x, c('hindfootLength', 'weight', 'logweight')])
siteregpoolsp_mam15_iucn <- lapply(reg_pool[sites15], function(x) mam_capture_sitemerge[x, c('taxonID')])

# Define entire continent as regional pool for each of the local and regional stats.
continentpool15 <- lapply(reg_pool[sites15], function(x) mam_capture_sitemerge[, c('hindfootLength', 'weight', 'logweight')])
continentsp15 <- lapply(reg_pool[sites15], function(x) mam_capture_sitemerge[, c('taxonID')])

# RUN THE OSTATS

task <- as.numeric(Sys.getenv('PBS_ARRAYID'))
print(task)

set.seed(27510)

if (task == 1)
{Ostats_bysite2015local <- Ostats(traits = traits_mam15, plots = factor(mam_capture_sitemerge$siteID), sp = factor(mam_capture_sitemerge$taxonID), reg_pool_traits = siteregpoollist_mam15_local, reg_pool_sp = siteregpoolsp_mam15_local, nperm = N_PERM, stat = STAT, shuffle_weights = SW)
save(Ostats_bysite2015local, file = '~/data/Ostats09Jan/Ostats_bysite2015_harmoniclocal_orignull.r')}

if (task == 2)
{Ostats_bysite2015iucn <- Ostats(traits = traits_mam15, plots = factor(mam_capture_sitemerge$siteID), sp = factor(mam_capture_sitemerge$taxonID), reg_pool_traits = siteregpoollist_mam15_iucn, reg_pool_sp = siteregpoolsp_mam15_iucn, nperm = N_PERM, stat = STAT, shuffle_weights = SW)
save(Ostats_bysite2015iucn, file = '~/data/Ostats09Jan/Ostats_bysite2015_harmoniciucn_orignull.r')}

if (task == 3)
{Ostats_bysite2015conti <- Ostats(traits = traits_mam15, plots = factor(mam_capture_sitemerge$siteID), sp = factor(mam_capture_sitemerge$taxonID), reg_pool_traits = continentpool15, reg_pool_sp = continentsp15, nperm = N_PERM, stat = STAT, shuffle_weights = SW)
save(Ostats_bysite2015conti, file = '~/data/Ostats09Jan/Ostats_bysite2015_harmonicconti_orignull.r')}
