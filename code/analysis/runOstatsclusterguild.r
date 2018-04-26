# Overlap Null Models (without T-stats) run on the mammal data. Run on cluster and split by guild.
# Author: QDR
# Project: NEON ITV
# Created: 11 Oct 2016 *from older script

# 25 Oct: Added code to get rid of the insectivorous mammals and the individuals not identified to species.
# 13 Oct: Add weighted median.
# 11 Oct: change to run remotely and split by guild. 
# Modified 10 Oct: Add regional O-stats calculation with continent as regional pool.

# Set number of null model iterations here and which stat is used.
STAT <- 'weighted.median'
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

# Here, get rid of the insectivores and the mammals not identified to species.
nottospecies <- c('DPSP','GLSP','MISP','NESP','PESP','PGSP','RESP')
mam_capture_sitemerge <- mam_capture_sitemerge %>% filter(Pineda_Main_food != 'Insectivore' & !(taxonID %in% nottospecies)) # Gets rid of ~700.

# Split into year+guild combination
i2014 <- mam_capture_sitemerge$year == 2014
i2015 <- mam_capture_sitemerge$year == 2015
iall <- rep(TRUE, length(i2014))

i2014gen <- mam_capture_sitemerge$year == 2014 & mam_capture_sitemerge$Pineda_Main_food == 'Generalist'
i2015gen <- mam_capture_sitemerge$year == 2015 & mam_capture_sitemerge$Pineda_Main_food == 'Generalist'
iallgen <- mam_capture_sitemerge$Pineda_Main_food == 'Generalist'

# Convert adult mammal data to trait-only format
traits_mam14 <- mam_capture_sitemerge[i2014, c('hindfootLength', 'weight', 'logweight')]
traits_mam15 <- mam_capture_sitemerge[i2015, c('hindfootLength', 'weight', 'logweight')]
traits_mamallyears <- mam_capture_sitemerge[, c('hindfootLength', 'weight', 'logweight')]

traits_gen14 <- mam_capture_sitemerge[i2014gen, c('hindfootLength', 'weight', 'logweight')]
traits_gen15 <- mam_capture_sitemerge[i2015gen, c('hindfootLength', 'weight', 'logweight')]
traits_genallyears <- mam_capture_sitemerge[iallgen, c('hindfootLength', 'weight', 'logweight')]

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

sites14gen <- unique(mam_capture_sitemerge$siteID[i2014gen])
sites15gen <- unique(mam_capture_sitemerge$siteID[i2015gen])
plots14gen <- unique(mam_capture_sitemerge$plotID[i2014gen])
plots15gen <- unique(mam_capture_sitemerge$plotID[i2015gen])
sitesallyearsgen <- unique(mam_capture_sitemerge$siteID[iallgen])
plotsallyearsgen <- unique(mam_capture_sitemerge$plotID[iallgen])

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

siteregpoollist_gen14_iucn <- lapply(reg_pool[sites14gen], function(x) mam_capture_sitemerge[x, c('hindfootLength', 'weight', 'logweight')])
siteregpoollist_gen15_iucn <- lapply(reg_pool[sites15gen], function(x) mam_capture_sitemerge[x, c('hindfootLength', 'weight', 'logweight')])
siteregpoollist_genallyears_iucn <- lapply(reg_pool[sitesallyearsgen], function(x) mam_capture_sitemerge[x, c('hindfootLength', 'weight', 'logweight')])

siteregpoolsp_gen14_iucn <- lapply(reg_pool[sites14gen], function(x) mam_capture_sitemerge[x, c('taxonID')])
siteregpoolsp_gen15_iucn <- lapply(reg_pool[sites15gen], function(x) mam_capture_sitemerge[x, c('taxonID')])
siteregpoolsp_genallyears_iucn <- lapply(reg_pool[sitesallyearsgen], function(x) mam_capture_sitemerge[x, c('taxonID')])


# Define entire continent as regional pool for each of the local and regional stats.
continentpool14 <- lapply(reg_pool[sites14], function(x) mam_capture_sitemerge[, c('hindfootLength', 'weight', 'logweight')])
continentpool15 <- lapply(reg_pool[sites15], function(x) mam_capture_sitemerge[, c('hindfootLength', 'weight', 'logweight')])
continentpoolallyears <- lapply(reg_pool, function(x) mam_capture_sitemerge[, c('hindfootLength', 'weight', 'logweight')])

continentsp14 <- lapply(reg_pool[sites14], function(x) mam_capture_sitemerge[, c('taxonID')])
continentsp15 <- lapply(reg_pool[sites15], function(x) mam_capture_sitemerge[, c('taxonID')])
continentspallyears <- lapply(reg_pool, function(x) mam_capture_sitemerge[, c('taxonID')])

continentpool14gen <- lapply(reg_pool[sites14gen], function(x) mam_capture_sitemerge[, c('hindfootLength', 'weight', 'logweight')])
continentpool15gen <- lapply(reg_pool[sites15gen], function(x) mam_capture_sitemerge[, c('hindfootLength', 'weight', 'logweight')])
continentpoolallyearsgen <- lapply(reg_pool[sitesallyearsgen], function(x) mam_capture_sitemerge[, c('hindfootLength', 'weight', 'logweight')])

continentsp14gen <- lapply(reg_pool[sites14gen], function(x) mam_capture_sitemerge[, c('taxonID')])
continentsp15gen <- lapply(reg_pool[sites15gen], function(x) mam_capture_sitemerge[, c('taxonID')])
continentspallyearsgen <- lapply(reg_pool[sitesallyearsgen], function(x) mam_capture_sitemerge[, c('taxonID')])

# RUN THE OSTATS

task <- as.numeric(Sys.getenv('PBS_ARRAYID'))
print(task)

if (task == 1)
{Ostats_bysite2014 <- Ostats(traits = traits_mam14, plots = factor(mam_capture_sitemerge$siteID[i2014]), sp = factor(mam_capture_sitemerge$taxonID[i2014]), reg_pool_traits = siteregpoollist_mam14_iucn, reg_pool_sp = siteregpoolsp_mam14_iucn, nperm = N_PERM, stat = STAT)
save(Ostats_bysite2014, file = '~/data/Ostats_bysite2014_wmed.r')}
if (task == 2)
{Ostats_bysite2015 <- Ostats(traits = traits_mam15, plots = factor(mam_capture_sitemerge$siteID[i2015]), sp = factor(mam_capture_sitemerge$taxonID[i2015]), reg_pool_traits = siteregpoollist_mam15_iucn, reg_pool_sp = siteregpoolsp_mam15_iucn, nperm = N_PERM, stat = STAT)
save(Ostats_bysite2015, file = '~/data/Ostats_bysite2015_wmed.r')}
if (task == 3)
{Ostats_bysiteallyears <- Ostats(traits = traits_mamallyears, plots = factor(mam_capture_sitemerge$siteID[iall]), sp = factor(mam_capture_sitemerge$taxonID[iall]), reg_pool_traits = siteregpoollist_mamallyears_iucn, reg_pool_sp = siteregpoolsp_mamallyears_iucn, nperm = N_PERM, stat = STAT)
save(Ostats_bysiteallyears, file = '~/data/Ostats_bysiteallyears_wmed.r')}

# Ostats_regional2014 <- Ostats_regional(traits = traits_mam14, plots = factor(mam_capture_sitemerge$siteID[i2014]), sp = factor(mam_capture_sitemerge$taxonID[i2014]), reg_pool_traits = siteregpoollist_mam14_iucn, reg_pool_sp = siteregpoolsp_mam14_iucn, nperm = N_PERM)
# Ostats_regional2015 <-  Ostats_regional(traits = traits_mam15, plots = factor(mam_capture_sitemerge$siteID[i2015]), sp = factor(mam_capture_sitemerge$taxonID[i2015]), reg_pool_traits = siteregpoollist_mam15_iucn, reg_pool_sp = siteregpoolsp_mam15_iucn, nperm = N_PERM)
# Ostats_regionalallyears <- Ostats_regional(traits = traits_mamallyears, plots = factor(mam_capture_sitemerge$siteID[iall]), sp = factor(mam_capture_sitemerge$taxonID[iall]), reg_pool_traits = siteregpoollist_mamallyears_iucn, reg_pool_sp = siteregpoolsp_mamallyears_iucn, nperm = N_PERM)

# Ostats_regcontinent2014 <- Ostats_regional(traits = traits_mam14, plots = factor(mam_capture_sitemerge$siteID[i2014]), sp = factor(mam_capture_sitemerge$taxonID[i2014]), reg_pool_traits = continentpool14, reg_pool_sp = continentsp14, nperm = N_PERM)
# Ostats_regcontinent2015 <-  Ostats_regional(traits = traits_mam15, plots = factor(mam_capture_sitemerge$siteID[i2015]), sp = factor(mam_capture_sitemerge$taxonID[i2015]), reg_pool_traits = continentpool15, reg_pool_sp = continentsp15, nperm = N_PERM)
# Ostats_regcontinentallyears <- Ostats_regional(traits = traits_mamallyears, plots = factor(mam_capture_sitemerge$siteID[iall]), sp = factor(mam_capture_sitemerge$taxonID[iall]), reg_pool_traits = continentpoolallyears, reg_pool_sp = continentspallyears, nperm = N_PERM)

#save(Ostats_regcontinent2014, Ostats_regcontinent2015, Ostats_regcontinentallyears, file = '~/data/Ostats_continentpool.r')



if (task==4) 
{genOstats_bysite2014 <- Ostats(traits = traits_gen14, plots = factor(mam_capture_sitemerge$siteID[i2014gen]), sp = factor(mam_capture_sitemerge$taxonID[i2014gen]), reg_pool_traits = siteregpoollist_gen14_iucn, reg_pool_sp = siteregpoolsp_gen14_iucn, nperm = N_PERM, stat = STAT)
save(genOstats_bysite2014, file='~/data/genOstats_bysite2014_wmed.r')}
if (task==5)
{genOstats_bysite2015 <- Ostats(traits = traits_gen15, plots = factor(mam_capture_sitemerge$siteID[i2015gen]), sp = factor(mam_capture_sitemerge$taxonID[i2015gen]), reg_pool_traits = siteregpoollist_mam15_iucn, reg_pool_sp = siteregpoolsp_mam15_iucn, nperm = N_PERM, stat = STAT)
save(genOstats_bysite2015, file='~/data/genOstats_bysite2015_wmed.r')}
if (task==6)
{genOstats_bysiteallyears <- Ostats(traits = traits_genallyears, plots = factor(mam_capture_sitemerge$siteID[iallgen]), sp = factor(mam_capture_sitemerge$taxonID[iallgen]), reg_pool_traits = siteregpoollist_mamallyears_iucn, reg_pool_sp = siteregpoolsp_mamallyears_iucn, nperm = N_PERM, stat = STAT)
save(genOstats_bysiteallyears, file='~/data/genOstats_bysiteallyears_wmed.r')}
if (task==7)
{genOstats_regional2014 <- Ostats_regional(traits = traits_gen14, plots = factor(mam_capture_sitemerge$siteID[i2014gen]), sp = factor(mam_capture_sitemerge$taxonID[i2014gen]), reg_pool_traits = siteregpoollist_mam14_iucn, reg_pool_sp = siteregpoolsp_mam14_iucn, nperm = N_PERM)
save(genOstats_regional2014, file='~/data/genOstats_regional2014.r')}
if (task==8)
{genOstats_regional2015 <-  Ostats_regional(traits = traits_gen15, plots = factor(mam_capture_sitemerge$siteID[i2015gen]), sp = factor(mam_capture_sitemerge$taxonID[i2015gen]), reg_pool_traits = siteregpoollist_mam15_iucn, reg_pool_sp = siteregpoolsp_mam15_iucn, nperm = N_PERM)
save(genOstats_regional2015, file='~/data/genOstats_regional2015.r')}
if (task==9)
{genOstats_regionalallyears <- Ostats_regional(traits = traits_genallyears, plots = factor(mam_capture_sitemerge$siteID[iallgen]), sp = factor(mam_capture_sitemerge$taxonID[iallgen]), reg_pool_traits = siteregpoollist_mamallyears_iucn, reg_pool_sp = siteregpoolsp_mamallyears_iucn, nperm = N_PERM)
save(genOstats_regionalallyears, file='~/data/genOstats_regionalallyears.r')}
if (task==10)
{genOstats_regcontinent2014 <- Ostats_regional(traits = traits_gen14, plots = factor(mam_capture_sitemerge$siteID[i2014gen]), sp = factor(mam_capture_sitemerge$taxonID[i2014gen]), reg_pool_traits = continentpool14gen, reg_pool_sp = continentsp14gen, nperm = N_PERM)
save(genOstats_regcontinent2014, file='~/data/genOstats_regcontinent2014.r')}
if (task==11)
{genOstats_regcontinent2015 <-  Ostats_regional(traits = traits_gen15, plots = factor(mam_capture_sitemerge$siteID[i2015gen]), sp = factor(mam_capture_sitemerge$taxonID[i2015gen]), reg_pool_traits = continentpool15gen, reg_pool_sp = continentsp15gen, nperm = N_PERM)
save(genOstats_regcontinent2015, file='~/data/genOstats_regcontinent2015.r')}
if (task==12)
{genOstats_regcontinentallyears <- Ostats_regional(traits = traits_genallyears, plots = factor(mam_capture_sitemerge$siteID[iallgen]), sp = factor(mam_capture_sitemerge$taxonID[iallgen]), reg_pool_traits = continentpoolallyearsgen, reg_pool_sp = continentspallyearsgen, nperm = N_PERM)
save(genOstats_regcontinentallyears, file='~/data/genOstats_regcontinentallyears.r')}

# no-ITV cases
if (task == 13)
{Ostats_bysite2014_noITV <- Ostats(traits = traits_mam14, plots = factor(mam_capture_sitemerge$siteID[i2014]), sp = factor(mam_capture_sitemerge$taxonID[i2014]), reg_pool_traits = siteregpoollist_mam14_iucn, reg_pool_sp = siteregpoolsp_mam14_iucn, nperm = N_PERM, stat = 'none')
save(Ostats_bysite2014_noITV, file = '~/data/Ostats_bysite2014_noitv.r')}
if (task == 14)
{Ostats_bysite2015_noITV <- Ostats(traits = traits_mam15, plots = factor(mam_capture_sitemerge$siteID[i2015]), sp = factor(mam_capture_sitemerge$taxonID[i2015]), reg_pool_traits = siteregpoollist_mam15_iucn, reg_pool_sp = siteregpoolsp_mam15_iucn, nperm = N_PERM, stat = 'none')
save(Ostats_bysite2015_noITV, file = '~/data/Ostats_bysite2015_noitv.r')}
if (task == 15)
{Ostats_bysiteallyears_noITV <- Ostats(traits = traits_mamallyears, plots = factor(mam_capture_sitemerge$siteID[iall]), sp = factor(mam_capture_sitemerge$taxonID[iall]), reg_pool_traits = siteregpoollist_mamallyears_iucn, reg_pool_sp = siteregpoolsp_mamallyears_iucn, nperm = N_PERM, stat = 'none')
save(Ostats_bysiteallyears_noITV, file = '~/data/Ostats_bysiteallyears_noitv.r')}