# Calculation of mammal phylogenetic diversity at NEON sites
# Author: QDR
# Project: NEON ITV
# Created: 22 Aug 2016
# Last modified: 12 Jan 2017

# Modified 20 Apr 2018: Specify paths

# Modified 12 Jan: Recalculate using new species list (only target species, just rodents)
# Modified 26 Aug: calculate PD by PLOT in addition to by site, also change to independent swap algorithm.


# Data loading and processing ---------------------------------------------

data_path <- '/mnt/research/neon'

# Load in the phylogeny for mammals

library(ape)
library(picante)
tol <- read.tree(file.path(data_path, 'external_data/raw_external_data/phylo/TimetreeOfLife2015.nwk'))

# Load mammal names
mammalTax <- read.csv('species_lists/mammalSpeciesList.csv', stringsAsFactors = FALSE)

# Load NEON data and calculate phylogenetic diversity indices
load(file.path(data_path, 'final_data/allorganismal_latest.r'))

# Process data to get communities for each year (same communities as are used for the Chao calculation)

library(dplyr)
library(lubridate)
library(reshape2)

# These are all the capture grids, with possibly multiple bouts per year
# Get rid of blank or "Other" taxonIDs, and two individuals with a possibly typo species code
# mam_comms <- mam_capture %>% mutate(gridID = substr(trapCoordinate,1,1), year = year(date)) %>%
#   filter(taxonID != '', taxonID != 'OTHE', taxonID != 'ROSP') %>%
#   group_by(year, siteID) %>%
#   do(comm = dcast(., formula = gridID ~ taxonID)[,-1, drop = FALSE]) # Get rid of date col and use n. rows as n. individuals
# 
# mam_comms_byyear <- do(mam_comms, abund = apply(.$comm, 2, sum))

# Make a single matrix for each year.
comm_mat <- function(x) {
  mat <- dcast(x, formula = siteID ~ taxonID)
  row.names(mat) <- mat$siteID
  mat[,-1, drop=FALSE]
}

# Matrix by plot not just site
comm_mat_byplot <- function(x) {
  mat <- dcast(x, formula = plotID ~ taxonID)
  row.names(mat) <- mat$plotID
  mat[,-1, drop=FALSE]
}

mam_use <- mam_capture %>% 
  mutate(year = year(date)) %>%
  filter(taxonProtocolCategory == 'target', taxonID != 'OTHE', year != 2016)
  
 # All years.
mam_use <- mam_capture %>% 
  mutate(year = year(date)) %>%
  filter(taxonProtocolCategory == 'target', taxonID != 'OTHE')
 

mam_mats <- mam_use %>%
  group_by(year) %>%
  do(comm = comm_mat(.))


# mam_mats <- mam_capture %>% mutate(year = year(date)) %>%
#   filter(taxonID != '', taxonID != 'OTHE', taxonID != 'ROSP') %>%
#   group_by(year) %>%
#   do(comm = comm_mat(.))

# To do by plot (Aug 26):
mam_mats <- mam_use %>% mutate(year = year(date)) %>%
  filter(taxonID != '', taxonID != 'OTHE', taxonID != 'ROSP') %>%
  group_by(year) %>%
  do(comm = comm_mat_byplot(.))

# Get the list of four-letter codes that are in the mammal matrices
allIDs <- Reduce(union, sapply(mam_mats$comm, names))
allIDs <- merge(data.frame(taxonID=allIDs), mammalTax[,c('taxonID','scientificName')], all.x=TRUE)
allIDs <- transform(allIDs, scientificName = sub(" ", "_", scientificName))

# Figure out which of the mammals in mammal community aren't in the tree, and add them to "tol" using phyTools
# It is a lot of unidentified species, and some Spermophilus species.

spp_not_in_tree <- allIDs$scientificName[!allIDs$scientificName %in% tol$tip.label]

# Extract the mammal clade from phylogeny
sum(mammalTax$scientificName %in% tol$tip.label) # Wrong format
sum(with(mammalTax, paste(genus, specificEpithet, sep='_')) %in% tol$tip.label) # 300 are in

#not_mammals <- which(!(tol$tip.label %in% allIDs$scientificName) & !grepl('Spermophilus', tol$tip.label))
not_mammals <- which(!(tol$tip.label %in% allIDs$scientificName))
mammaltol <- drop.tip(tol, tip=not_mammals) 

# Add tips for individuals that are only identified to genus, and add Spermophilus.
# Edit: Spermophilus isn't needed anymore. Not target species.
library(phytools)
for (i in spp_not_in_tree) mammaltol <- add.species.to.genus(tree = mammaltol, species = i, where = 'root')
#for (i in spp_not_in_tree[10:13]) mammaltol <- add.species.to.genus(tree = mammaltol, species = i, where = 'random')

# Prune out the other spermophilus entries
mammaltol <- drop.tip(mammaltol, tip=which(!(mammaltol$tip.label %in% allIDs$scientificName)))

# Change the tip labels of the mammal phylogeny to the taxonIDs
mammaltol$tip.label <- as.character(allIDs$taxonID)[match(mammaltol$tip.label, allIDs$scientificName)]

# Save phylogeny for later use running phylogenetic corrections or other analyses.
save(mammaltol, file = 'mammaltol.r')
write.tree(mammaltol, 'mammaltol.nwk')

# Calculate pd metrics ----------------------------------------------------

# PD
mc2014 <- mam_mats$comm[[2]]
mc2015 <- mam_mats$comm[[3]]

tree2014 <- prune.sample(mc2014, mammaltol)
tree2015 <- prune.sample(mc2015, mammaltol)

pd2014 <- pd(mc2014, tree2014)
pd2015 <- pd(mc2015, tree2015)

# MPD and MNTD
dist2014 <- cophenetic(tree2014)
dist2015 <- cophenetic(tree2015)

mpd2014 <- ses.mpd(mc2014, dist2014, null.model = 'independentswap', abundance.weighted = TRUE, runs = 999, iterations = 1000)
mpd2015 <- ses.mpd(mc2015, dist2015, null.model = 'independentswap', abundance.weighted = TRUE, runs = 999, iterations = 1000)

mntd2014 <- ses.mntd(mc2014, dist2014, null.model = 'independentswap', abundance.weighted = TRUE, runs = 999, iterations = 1000)
mntd2015 <- ses.mntd(mc2015, dist2015, null.model = 'independentswap', abundance.weighted = TRUE, runs = 999, iterations = 1000)

names(pd2014) <- paste(names(pd2014), '2014', sep='')
names(pd2015) <- paste(names(pd2015), '2015', sep='')
names(mpd2014) <- paste(names(mpd2014), '2014', sep='')
names(mpd2015) <- paste(names(mpd2015), '2015', sep='')
names(mntd2014) <- paste(names(mntd2014), '2014', sep='')
names(mntd2015) <- paste(names(mntd2015), '2015', sep='')

pd2014$siteID <- row.names(pd2014)
pd2015$siteID <- row.names(pd2015)
mpd2014$siteID <- row.names(mpd2014)
mpd2015$siteID <- row.names(mpd2015)
mntd2014$siteID <- row.names(mntd2014)
mntd2015$siteID <- row.names(mntd2015)

# To do by plot (Aug 26)
pd2014$plotID <- row.names(pd2014)
pd2015$plotID <- row.names(pd2015)
mpd2014$plotID <- row.names(mpd2014)
mpd2015$plotID <- row.names(mpd2015)
mntd2014$plotID <- row.names(mntd2014)
mntd2015$plotID <- row.names(mntd2015)

mammalPD <- Reduce(full_join, list(pd2014, pd2015, mpd2014, mpd2015, mntd2014, mntd2015))
save(mammalPD, file = 'mammalPDobject.r')
# if done by plot
# This wasn't redone in January. Only the by-site was done.
save(mammalPD, file = 'mammalPDbyplotobject.r')

# Comparison of phylogenetic diversity metrics ----------------------------

# Compare null models within MPD and within MNTD

mpds <- list()
mntds <- list()
nms <- c("taxa.labels", "richness", "frequency", "sample.pool", "phylogeny.pool", "independentswap", "trialswap")

for (nm in nms) {
  mpds <- within(mpds, assign(nm, ses.mpd(mc2015, dist2015, null.model = nm, abundance.weighted = TRUE, runs = 999)$mpd.obs.z))
  mntds <- within(mntds, assign(nm, ses.mntd(mc2015, dist2015, null.model = nm, abundance.weighted = TRUE, runs = 999)$mntd.obs.z))
}

library(GGally)
ggpairs(do.call('cbind', mpds))
ggpairs(do.call('cbind', mntds))

# Compare the same null model algorithm across MPD and MNTD
qplot(mpds$taxa.labels, mntds$taxa.labels)
qplot(mpds$independentswap, mntds$independentswap)
