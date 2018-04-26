# Mammal analysis: Split by phylogenetic group (family)
# Author: QDR
# Project: NEON ITV
# Created: 22 July 2016
# Last modified: 20 Apr 2018

# Modified 20 Apr: specify path to data
data_path <- '/mnt/research/neon'

load('allorganismal2016Jun29.r') # New organismal data 
load('allsiteplot2016Jun21.r') # Site covariates

# Load full taxonomic categories
mammalTax <- read.csv('species_lists/mammalSpeciesList.csv', stringsAsFactors = FALSE)

# Load mammal phylogeny
library(ape)
timetree <- read.tree(file.path(data_path, 'external_data/raw_external_data/phylo/TimetreeOfLife2015.nwk'))

# Check where they match
mammalmatches <- with(mammalTax, paste(genus, specificEpithet, sep = '_')) %in% timetree$tip.label
mammalTax$scientificName[mammalmatches]

# Create only mammal phylogeny
to_delete <- timetree$tip.label[!timetree$tip.label %in% with(mammalTax, paste(genus, specificEpithet, sep = '_'))]
mammaltree <- drop.tip(phy=timetree, tip=to_delete)

# Further subset the tree to only the mammal species that actually appear in the surveys.
unique(mam_capture$scientificName) %in% mammalTax$scientificName[mammalmatches] # Almost all area present, except unknown species and the genus Spermophilus (family Sciuridae, tribe Marmotini)
unique(mam_capture$scientificName)[!unique(mam_capture$scientificName) %in% mammalTax$scientificName[mammalmatches]]

# Spermophilus should go in with Sciuridae

# Basic analysis: Sort by family.
mam_capture <- merge(mam_capture, mammalTax[,c('taxonID', 'family')], sort=FALSE, all.x=TRUE, all.y=FALSE)
