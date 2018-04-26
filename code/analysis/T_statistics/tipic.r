# Calculations of TIP/TIC statistics for different NEON taxa
# Author: QDR
# Project: NEON ITV
# Created: 1 Aug 2016
# Last modified: 21 Sep 2016

# Modified 21 Sep: Add log-transformation of body weights to the t-stat by subplot calculation, get rid of extraneous spp, add juveniles back into the analysis, and split analysis by guild too.
# Modified 24 Aug: Use fixed mammal plotIDs to run t-stats within sites!
# Modified 23 Aug: Attempt to classify by mammal grid so that we can look at within-site variation. Not much luck yet.
# Modified 19 Aug: Fix a problem with the indexing of regional species pool so that null models now work correctly.
# Modified 18 Aug: Use IUCN range maps to find regional species pool for each site.
# Modified 15 Aug: Split analysis into separate T-statistics for 2014 and 2015.
# Modified 11 Aug, 12 Aug: add the newest mammal data.

library(cati) # Package for calculating Tstats


# example data from package -----------------------------------------------



# Default call of Tstats
# Tstats(traits, ind.plot, sp, SE = 0, reg.pool = NULL, SE.reg.pool = NULL, nperm = 99, printprogress = TRUE)
# traits: matrix, individuals x traits
# ind.plot: vector, factor with plots
# sp: vector, factor with species
# SE: standard errors of each trait
# reg.pool: If there are more individuals with traits in the regional pool, we can use these

# Example
data(finch.ind)
# Loads traits.finch, ind.plot.finch, and sp.finch
res.finch <- Tstats(traits=traits.finch, ind.plot=ind.plot.finch, sp=sp.finch, nperm=9)

str(res.finch)
res.finch$Tstats$T_IP.IC
plot(res.finch)
plot(res.finch, type = "simple")
plot(res.finch, type = "simple_range")
plot(res.finch, type = "barplot")
plot(res.finch, type = "bysites")
plot(res.finch, type = "bytraits")

sum_Tstats(res.finch)
barplot(res.finch)


# Custom theme (from rasterVis package)
require(rasterVis)

my.theme <- BuRdTheme()
# Customize the colorkey
my.ckey <- list(col = my.theme$regions$col)

# Heat map of standardized effect sizes of TIPIC, as compared to the null model, on a site x trait matrix.

levelplot(t(ses(res.finch$Tstats$T_IP.IC,res.finch$Tstats$T_IP.IC_nm)$ses), 
          colorkey = my.ckey, par.settings = my.theme,border = "black")

# Create a fake regional species pool using some random samples from individuals in the Galapagos finch trait matrix.

reg.p <- rbind(traits.finch, traits.finch[sample(1:2000,300), ])

res.finch2 <- Tstats(traits.finch, ind.plot = ind.plot.finch, 
                     sp = sp.finch, reg.pool=reg.p, nperm = 9)	

plot(as.listofindex(list(res.finch,res.finch2)))

# Create a situation in which each island has a different regional species pool. Pass this list in to the reg.pool argument, instead of just a matrix as in the example above.


#### Use a different regional pool for each communities
#create a random regional pool for each communities for the example
list.reg.p <- list(
  traits.finch[sample(1:290,200), ], traits.finch[sample(100:1200,300), ], 
  traits.finch[sample(100:1500, 1000), ], traits.finch[sample(300:800,300), ],
  traits.finch[sample(1000:2000, 500), ], traits.finch[sample(100:900, 700), ] )

# Warning: the regional pool need to be larger than the observed communities
table(ind.plot.finch)
# For exemple, the third community need a regional pool of more than 981 individuals

res.finch3 <- Tstats(traits.finch, ind.plot = ind.plot.finch, 
                     sp = sp.finch, reg.pool=list.reg.p, nperm = 9, print = FALSE)	

plot(as.listofindex(list(res.finch, res.finch2, res.finch3)))	

# Show differences when accounting, and not accounting, for measurement error (SE of measurement set to 5 to model measurement error)

#### Use the standard errors of measure in the analysis (argument SE)
## Not run: 
res.finch.SE0 <- Tstats(traits.finch, ind.plot = ind.plot.finch, 
                        sp = sp.finch, SE = 0, print = FALSE)

res.finch.SE5 <- Tstats(traits.finch, ind.plot = ind.plot.finch, 
                        sp = sp.finch, SE = 5, print = FALSE)

plot(as.listofindex(list(res.finch.SE0, res.finch.SE5)))


# Calculate indices on mammal data ------------------------------------------

# Load neon data
load('allorganismal_latest.r') # Organismal data
load('allsiteplot_latest.r') # Site covariates

# Correct dates in mammal capture data frame. No longer necessary, I think
#date1 <- ymd(mam_capture$date)
#date1[is.na(date1)] <- mdy(mam_capture$date[is.na(date1)])
#mam_capture$date <- date1

# Here, subset out weasels, hares, and shrews.
# Get Rodentia only.

mammalTax <- read.csv('C:/Users/Q/Dropbox/neon/protocols/taxonomy/NEON_mam_taxonomy.csv')
mammalTraits <- read.csv('C:/Users/Q/Dropbox/neon/data/external_datasets/NEON_miscmammaltraits.csv')

# Add column to code for feeding guild, and run separate analysis for each guild.


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

mam_capture <- left_join(mam_capture, mammalGuilds[-c(89,118,62,63),], by = 'scientificName')

mam_noID <- mam_capture[mam_capture$individualandtag=='', c('year', 'siteID','taxonID','family', 'individualandtag','hindfootLength','earLength','tailLength','totalLength','weight','sex','lifeStage','Pineda_Main_food')]

mam_grp <- filter(mam_capture, individualandtag != '') %>%
  group_by(year, siteID, taxonID, family, individualandtag)

mam_byindiv <- summarize_at(mam_grp, vars(hindfootLength, earLength, tailLength, totalLength, weight), median, na.rm=TRUE)
mam_byindiv_class <- do(mam_grp, data.frame(sex = .$sex[1], lifeStage = names(sort(table(.$lifeStage),decreasing=TRUE))[1], Pineda_Main_food=.$Pineda_Main_food[1]))

mam_byindiv <- cbind(as.data.frame(mam_byindiv), as.data.frame(mam_byindiv_class[,c('sex', 'lifeStage','Pineda_Main_food')]))                  
mam_byindiv <-rbind(mam_byindiv, mam_noID)


# Merge covariates with mammal data for plotting
#mam_capture_merge <- merge(mam_byindiv, neonplotdata[,-c(1:6,9:10)], sort=FALSE, all.x=TRUE) # Does not work.
mam_capture_sitemerge <- merge(mam_byindiv, neonsitedata[,-c(3:5)], sort=FALSE, all.x=TRUE)
#mam_adults <- subset(mam_capture_sitemerge, lifeStage=='adult')

# Use juveniles as well as adults! (21 Sep)
mam_adults <- mam_capture_sitemerge




# Split into years
i2014 <- mam_adults$year == 2014
i2015 <- mam_adults$year == 2015

# Convert adult mammal data to cati format
traits_mam14 <- mam_adults[i2014, c('hindfootLength', 'earLength', 'tailLength', 'totalLength', 'weight')]
traits_mam15 <- mam_adults[i2015, c('hindfootLength', 'earLength', 'tailLength', 'totalLength', 'weight')]
# ind.plot is just the siteID, and sp is the taxonID

# Calculate T-statistics for the mammals

# This uses the entire country as the regional species pool for each site.
tstats_mam14 <- Tstats(traits = traits_mam14, ind.plot = factor(mam_adults$siteID[i2014]), sp = factor(mam_adults$taxonID[i2014]), nperm = 999)
tstats_mam15 <- Tstats(traits = traits_mam15, ind.plot = factor(mam_adults$siteID[i2015]), sp = factor(mam_adults$taxonID[i2015]), nperm = 999)

# Get a regional species pool for each site, as defined by the "domain" that it's in.
domaincombs <- unique(neonplotdata[,c('domainID','siteID')])
mam_adults <- merge(mam_adults, domaincombs, all.x=TRUE, sort=FALSE) # With 2015 data, now most domains have 2 sites, some have 3

# The traits in the regional species pool incorporate all the years' data.
mam_domains <- split(mam_adults, mam_adults$domainID)
domains14 <- unique(mam_adults$domainID[i2014])
domains15 <- unique(mam_adults$domainID[i2015])
sites14 <- unique(mam_adults$siteID[i2014])
sites15 <- unique(mam_adults$siteID[i2015])


regpoollist_mam14 <- regpoollist_mam15 <- list()
for (i in levels(factor(mam_adults$siteID))) {
  if (i %in% sites14) regpoollist_mam14[[length(regpoollist_mam14)+1]] <- mam_domains[[domaincombs$domainID[domaincombs$siteID==i]]]
  if (i %in% sites15) regpoollist_mam15[[length(regpoollist_mam15)+1]] <- mam_domains[[domaincombs$domainID[domaincombs$siteID==i]]]
}

tstats_mam_regpools14 <- Tstats(traits = traits_mam14, ind.plot = factor(mam_adults$siteID[i2014]), sp = factor(mam_adults$taxonID[i2014]), reg.pool = regpoollist_mam14, nperm = 999)
tstats_mam_regpools15 <- Tstats(traits = traits_mam15, ind.plot = factor(mam_adults$siteID[i2015]), sp = factor(mam_adults$taxonID[i2015]), reg.pool = regpoollist_mam15, nperm = 999)

# Improved regional species pool with the iucn range maps

iucn <- read.csv('C:/Users/Q/Dropbox/neon/data/external_datasets/IUCN_mammal_ranges.csv')
mamnames <- unique(mam_capture$scientificName) 
mamnames[mamnames %in% iucn$binomial]
mamnames[!mamnames %in% iucn$binomial] # These include taxa not identified to species as well as two additional species: Tamias alpinus and Rattus rattus. Tamias alpinus should only be in California sites (although it is recorded as a single individual from Utah), and Rattus rattus should be in all sites. Put all the genus-level traits into all the regional species pools.

# Fortunately, only 470 rows in the mam_capture dataframe are not accounted for in IUCN.
sum(mam_capture$scientificName %in% mamnames[!mamnames %in% iucn$binomial])

nomatchspp <- mamnames[!mamnames %in% iucn$binomial]
nomatchspp <- nomatchspp[c(1:5,7:11,13)] # Only the ones from all sites

# Include all traits from all individuals whose species id matches species at the site, plus all traits from the individuals whose species id matches species with ranges overlapping the site (in case IUCN ranges are off)

reg_pool <- list()

for (site in unique(mam_capture$siteID)) {
  iucn_spp <- as.character(iucn$binomial[iucn[,site]])
  neon_spp <- mam_capture$scientificName[mam_capture$siteID == site]
  regional_spp <- unique(c(iucn_spp, neon_spp, nomatchspp))
  regional_spp_id <- unique(mam_capture$taxonID[mam_capture$scientificName %in% regional_spp])
  reg_pool <- within(reg_pool, assign(site, value=mam_adults$taxonID %in% regional_spp_id))
}

# This fix corrects a problem with the null models (adding column name selection)
regpoollist_mam14_iucn <- lapply(reg_pool[sites14], function(x) mam_adults[x, c('hindfootLength', 'earLength', 'tailLength', 'totalLength', 'weight')])
regpoollist_mam15_iucn <- lapply(reg_pool[sites15], function(x) mam_adults[x, c('hindfootLength', 'earLength', 'tailLength', 'totalLength', 'weight')])

tstats_mam_regpools14iucn <- Tstats(traits = traits_mam14, ind.plot = factor(mam_adults$siteID[i2014]), sp = factor(mam_adults$taxonID[i2014]), reg.pool = regpoollist_mam14_iucn, nperm = 999)
tstats_mam_regpools15iucn <- Tstats(traits = traits_mam15, ind.plot = factor(mam_adults$siteID[i2015]), sp = factor(mam_adults$taxonID[i2015]), reg.pool = regpoollist_mam15_iucn, nperm = 999)

# Save and load tstats objects as they take a very long time to generate
save(tstats_mam14, tstats_mam15, tstats_mam_regpools14, tstats_mam_regpools15, tstats_mam_regpools14iucn, tstats_mam_regpools15iucn, file='C:/Users/Q/Dropbox/neon/code/tstatsmam.r')

# Just regpools ojbect (15 aug)
#save(tstats_mam_regpools, file='C:/Users/Q/Dropbox/neon/code/tstatsmamregpools.r')


# Run Tstats by subplot (for mixed model) ---------------------------------

# 0. Process by plotID (mostly the same as above except plotID is added as a grouping variable.)

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
#mam_capture_merge <- merge(mam_byindiv, neonplotdata[,-c(1:6,9:10)], sort=FALSE, all.x=TRUE) # Does not work.
mam_capture_sitemerge <- merge(mam_byindiv, neonsitedata[,-c(3:5)], sort=FALSE, all.x=TRUE)
mam_adults <- subset(mam_capture_sitemerge, lifeStage=='adult') %>% mutate(logweight = log(weight)) # Add log transformation

# Split into years
i2014 <- mam_adults$year == 2014
i2015 <- mam_adults$year == 2015

# Convert adult mammal data to cati format
traits_mam14 <- mam_adults[i2014, c('hindfootLength', 'earLength', 'tailLength', 'totalLength', 'weight', 'logweight')]
traits_mam15 <- mam_adults[i2015, c('hindfootLength', 'earLength', 'tailLength', 'totalLength', 'weight', 'logweight')]

# 1. Define IUCN Regional Pools

iucn <- read.csv('C:/Users/Q/Dropbox/neon/data/external_datasets/IUCN_mammal_ranges.csv')
mamnames <- unique(mam_capture$scientificName) 
mamnames[mamnames %in% iucn$binomial]
mamnames[!mamnames %in% iucn$binomial] # These include taxa not identified to species as well as two additional species: Tamias alpinus and Rattus rattus. Tamias alpinus should only be in California sites (although it is recorded as a single individual from Utah), and Rattus rattus should be in all sites. Put all the genus-level traits into all the regional species pools.

# Fortunately, only 470 rows in the mam_capture dataframe are not accounted for in IUCN.
sum(mam_capture$scientificName %in% mamnames[!mamnames %in% iucn$binomial])

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

plotregpoollist_mam14_iucn <- lapply(reg_pool[plots14], function(x) mam_adults[x, c('hindfootLength', 'earLength', 'tailLength', 'totalLength', 'weight', 'logweight')])
plotregpoollist_mam15_iucn <- lapply(reg_pool[plots15], function(x) mam_adults[x, c('hindfootLength', 'earLength', 'tailLength', 'totalLength', 'weight', 'logweight')])

# 2. Define site-only Regional Pools

plotregpoollist_mam14_siteonly <- lapply(plots14, function(x) {
  spp <- unique(mam_adults$taxonID[mam_adults$siteID == substr(x, 1, 4)])
  mam_adults[mam_adults$taxonID %in% spp, c('hindfootLength', 'earLength', 'tailLength', 'totalLength', 'weight', 'logweight')]
})
plotregpoollist_mam15_siteonly <- lapply(plots15, function(x) {
  spp <- unique(mam_adults$taxonID[mam_adults$siteID == substr(x, 1, 4)])
  mam_adults[mam_adults$taxonID %in% spp, c('hindfootLength', 'earLength', 'tailLength', 'totalLength', 'weight', 'logweight')]
})

names(plotregpoollist_mam14_siteonly) <- plots14
names(plotregpoollist_mam15_siteonly) <- plots15

# 2b. Modification Aug 25: try excluding Healy since it has a very small number of individuals
# This only is relevant for the 2015 site-only regional pool.
# Only 16 rows are there.
i2015b <- mam_adults$year == 2015 & mam_adults$siteID != "HEAL"
traits_mam15b <- mam_adults[i2015b, c('hindfootLength', 'earLength', 'tailLength', 'totalLength', 'weight', 'logweight')]
sites15b <- unique(mam_adults$siteID[i2015b])
plots15b <- unique(mam_adults$plotID[i2015b])

plotregpoollist_mam15_siteonly <- lapply(plots15b, function(x) {
  spp <- unique(mam_adults$taxonID[mam_adults$siteID == substr(x, 1, 4)])
  mam_adults[mam_adults$taxonID %in% spp, c('hindfootLength', 'earLength', 'tailLength', 'totalLength', 'weight', 'logweight')]
})

# 3. Calculate T-stats for each year and each regional pool type

tstats_mam_plots14_iucnpool <- Tstats(traits = traits_mam14, ind.plot = factor(mam_adults$plotID[i2014]), sp = factor(mam_adults$taxonID[i2014]), reg.pool = plotregpoollist_mam14_iucn, nperm = 999)
tstats_mam_plots15_iucnpool <- Tstats(traits = traits_mam15, ind.plot = factor(mam_adults$plotID[i2015]), sp = factor(mam_adults$taxonID[i2015]), reg.pool = plotregpoollist_mam15_iucn, nperm = 999)
tstats_mam_plots14_sitepool <- Tstats(traits = traits_mam14, ind.plot = factor(mam_adults$plotID[i2014]), sp = factor(mam_adults$taxonID[i2014]), reg.pool = plotregpoollist_mam14_siteonly, nperm = 999)
#tstats_mam_plots15_sitepool <- Tstats(traits = traits_mam15, ind.plot = factor(mam_adults$plotID[i2015]), sp = factor(mam_adults$taxonID[i2015]), reg.pool = plotregpoollist_mam15_siteonly, nperm = 999)
tstats_mam_plots15_sitepool <- Tstats(traits = traits_mam15b, ind.plot = factor(mam_adults$plotID[i2015b]), sp = factor(mam_adults$taxonID[i2015b]), reg.pool = plotregpoollist_mam15_siteonly, nperm = 999)

save(tstats_mam_plots14_iucnpool, tstats_mam_plots15_iucnpool, tstats_mam_plots14_sitepool, tstats_mam_plots15_sitepool, file = 'C:/Users/Q/Dropbox/neon/code/tstatsmambyplot.r')

# Plot mammal Tstats ------------------------------------------------------



load('C:/Users/Q/Dropbox/neon/code/tstatsmam.r')

# Print some data visualizations.
plot(tstats_mam_regpools)

tstats_mam_regpools$Tstats$T_IP.IC
plot(tstats_mam_regpools, type = "simple")
plot(tstats_mam_regpools, type = "simple_range")
plot(tstats_mam_regpools, type = "barplot")
plot(tstats_mam_regpools, type = "bysites")
plot(tstats_mam_regpools, type = "bytraits")

sum_Tstats(tstats_mam_regpools)
barplot(tstats_mam_regpools)


# Custom theme (from rasterVis package)
require(rasterVis)

my.theme <- BuRdTheme()
# Customize the colorkey
my.ckey <- list(col = my.theme$regions$col)

# Heat map of standardized effect sizes of TIPIC, as compared to the null model, on a site x trait matrix.

levelplot(t(ses(tstats_mam_regpools$Tstats$T_IP.IC,tstats_mam_regpools$Tstats$T_IP.IC_nm)$ses), 
          colorkey = my.ckey, par.settings = my.theme,border = "black")

# Plot both types of t-stat against each other
plot(as.listofindex(tstats_mam, tstats_mam_regpools))


# Generate phenology data -------------------------------------------

phen_fg <- read.csv('C:/Users/Q/Dropbox/neon/data/pheno_plants_funcgrp.csv')
phe_status <- merge(phe_status, phen_fg[,-1], sort=FALSE, all.x=TRUE, all.y=FALSE)

# Convert phenophase names to comparable numbers.

phase_key <- c('Breaking leaf buds'=1, 'Colored leaves'=5, 'Emerging needles'=1, 'Falling leaves'=6, 'Increasing leaf size'=2, 'Initial growth'=1, 'Leaves'=3, 'Open flowers'=4, 'Open pollen cones'=4, 'Young leaves'=2, 'Young needles'=2)

phe_status$Key <- phase_key[phe_status$phenophaseName]

# Convert phenophase percentage classes to the midpoints of those classes
intensity_percentages <- c('< 5%'=0.025,
                           '5-24%'=0.15,
                           '<25%'=0.125,
                           '25-49%'=0.375,
                           '50-74%'=0.625,
                           '75-94%'=0.85,
                           '>= 95%'=0.975)

# Convert phenophase raw number classes to the midpoints of those classes
intensity_numbers <- c('< 3'=1,
                       '3 to 10'=6.5,
                       '11 to 100'=55.5,
                       '101 to 1000'=550.5,
                       '1001 to 10000'=5500.5)

# Scale this also to a zero-one scale by dividing by its max.
intensity_numbers <- intensity_numbers/max(intensity_numbers)

# Convert the categorical values to numerical.
pctvalues <- intensity_percentages[phe_status$phenophaseIntensity]
numvalues <- intensity_numbers[phe_status$phenophaseIntensity]

# Make a new column in phe_status that includes all the numerical values.
phe_status$intensityValue <- pmin(pctvalues, numvalues, na.rm=TRUE)

# Find the onset, end, and length for the phenophases for each individual (is it possible to do this for intensity)?

library(lubridate)
library(plyr)

phe_df <- phe_status[,c('domainID','siteID','plotID','date','dayOfYear','individualID','taxonID','func_grp','phenophaseStatus','Key','intensityValue')]

phenotable <- ddply(phe_df, .(siteID, individualID, taxonID, func_grp, year(date)), function(x) {
  onset <- rep(NA, 6)
  final <- rep(NA, 6)
  
  for (k in 1:6) {
    days_k <- x$dayOfYear[x$Key==k & x$phenophaseStatus=='yes']
    onset[k] <- ifelse(length(days_k>0), min(days_k), NA)
    final[k] <- ifelse(length(days_k>0), max(days_k), NA)
  }
  names(onset) <- paste('onset', 1:6, sep='_')
  names(final) <- paste('end', 1:6, sep='_')
  
  return(c(onset, final))
})

library(reshape2)
phenotable_onsets <- melt(phenotable[,1:11], id.vars=1:5, value.name='day_onset')
phenotable_ends <- melt(phenotable[,c(1:5, 12:17)], id.vars=1:5, value.name='day_end')
phenotable_onsets$phase <- sapply(strsplit(as.character(phenotable_onsets$variable), '_'), '[', 2)
phenotable_ends$phase <- sapply(strsplit(as.character(phenotable_ends$variable), '_'), '[', 2)

phenotable_long <- cbind(phenotable_onsets, day_end=phenotable_ends$day_end)

# Get a single number to be the individual ID for each species.
indivnos <- ddply(phenotable_long, .(siteID), function(x) data.frame(siteID=x$siteID, individualID=x$individualID, indivno=as.numeric(factor(x$individualID))))
phenotable_long$indivno <- indivnos$indivno[match(phenotable_long$individualID, indivnos$individualID)]

phenotable_long$phase_length <- with(phenotable_long, day_end - day_onset)
names(phenotable_long)[5] <- 'year'

# Calculate variability in different aspects of phenology
library(plyr)
phenotable_sd <- ddply(phenotable_long, .(siteID, taxonID, func_grp, year, phase), colwise(sd, .(day_onset, day_end, phase_length), na.rm=TRUE))
phenotable_n <- ddply(phenotable_long, .(siteID, taxonID, func_grp, year, phase), colwise(function(x) sum(!is.na(x)), .(day_onset, day_end, phase_length)))
phenotable_sd <- transform(phenotable_sd, 
                           n_onset = phenotable_n$day_onset,
                           n_end = phenotable_n$day_end,
                           n_phase_length = phenotable_n$phase_length
)

# add environmental covariates to phenotable_sd
phenotable_sd <- merge(phenotable_sd, neonsitedata, sort=FALSE, all.x=TRUE, all.y=FALSE)

# Indices with phenology data ---------------------------------------------

# Use phenotable (non-melted)

# Get only valid records
phenotable_valid <- subset(phenotable, !is.na(taxonID) & taxonID!='')

# Figure out where there are replicated individuals.
phenotable_valid <- subset(phenotable_valid, !(is.na(onset_1) & is.na(onset_2) & is.na(onset_3) & is.na(onset_4) & is.na(onset_5) & is.na(onset_6)))
nindivs <- table(phenotable_valid$individualID)

phenotable_valid[phenotable_valid$individualID %in% names(nindivs[nindivs>1]), ]

# Get rid of the duplicated ones that do not have a budbreak date, which should take care of the rest.
toremove <- row.names(phenotable_valid[phenotable_valid$individualID %in% names(nindivs[nindivs>1]) & is.na(phenotable_valid$onset_1), ])
phenotable_valid <- phenotable_valid[!(row.names(phenotable_valid) %in% toremove),]

# Use only the onset dates for the phenological traits
# Put in cati format
traits_phe <- phenotable_valid[,c('onset_1','onset_2','onset_3','onset_4','onset_5','onset_6')]

# Get a regional species pool for each site, as defined by the "domain" that it's in.
domaincombs <- unique(neonplotdata[,c('domainID','siteID')])
phenotable_valid <- merge(phenotable_valid, domaincombs, all.x=TRUE, sort=FALSE)

phe_domains <- split(traits_phe, phenotable_valid$domainID)
regpoollist_phe <- list()
for (i in levels(factor(phenotable_valid$siteID))) regpoollist_phe[[length(regpoollist_phe)+1]] <- phe_domains[[domaincombs$domainID[domaincombs$siteID==i]]]

# Calculate T-statistics with the entire country as a regional species pool.
tstats_phe <- Tstats(traits = traits_phe, ind.plot = factor(phenotable_valid$siteID), sp = factor(phenotable_valid$taxonID), nperm = 999)

tstats_phe_regpools <- Tstats(traits = traits_phe, ind.plot = factor(phenotable_valid$siteID), sp = factor(phenotable_valid$taxonID), reg.pool = regpoollist_phe, nperm = 999)

# Save and load tstats objects as they take a very long time to generate
save(tstats_phe, tstats_phe_regpools, file='C:/Users/Q/Dropbox/neon/code/tstatsphe.r')

load('C:/Users/Q/Dropbox/neon/code/tstatsphe.r')

# Visualize
plot(tstats_phe_regpools)

tstats_phe_regpools$Tstats$T_IP.IC
plot(tstats_phe_regpools, type = "simple")
plot(tstats_phe_regpools, type = "simple_range")
plot(tstats_phe_regpools, type = "barplot")
plot(tstats_phe_regpools, type = "bysites")
plot(tstats_phe_regpools, type = "bytraits")

sum_Tstats(tstats_phe_regpools)
barplot(tstats_phe_regpools)


# Custom theme (from rasterVis package)
require(rasterVis)

my.theme <- BuRdTheme()
# Customize the colorkey
my.ckey <- list(col = my.theme$regions$col)

# Heat map of standardized effect sizes of TIPIC, as compared to the null model, on a site x trait matrix.

levelplot(t(ses(tstats_phe_regpools$Tstats$T_IP.IC,tstats_phe_regpools$Tstats$T_IP.IC_nm)$ses), 
          colorkey = my.ckey, par.settings = my.theme,border = "black")
