# Calculation of functional diversity of mammals.
# Author: QDR
# Project: NEON ITV
# Created: 29 June 2016
# Last modified: 21 July 2016

data_path <- '/mnt/research/neon'

# Modified on 30 June: new mammal data was added to NEON database

# Load data.
#load('allorganismal2016Feb25.r') # Organismal data
load(file.path(data_path, 'final_data/allorganismal2016Jun29.r')) # New organismal data
load(file.path(data_path, 'allsiteplot2016Jun21.r')) # Site covariates

# First step is to determine what constitutes an assemblage. Each site has multiple traps set over multiple nights.
# We can split it by site, plot, and night (100 traps)

library(plyr)
richness_bynight <- ddply(mam_capture, .(siteID, plotID, date), summarize,
                          richness = length(unique(taxonID)),
                          abundance = length(taxonID))

# Need to do rarefaction because of the greatly differing sampling effort, and number of individuals caught, across sites.

# Rarefaction Curve for plots at each site

# As of 21 July, this code is broken because the plotIDs in the mammal data are all set to NULL.

# Create site x species matrix for sites/sitexplots/sitexplotxnights

# Remove DSNY because it only has 6 individuals, respectively. (don't need to do this any more)
#mam_subset <- subset(mam_capture, !siteID %in% c('DSNY'))
mam_subset <- mam_capture

M_sites <- with(mam_subset, table(siteID, taxonID))
M_siteplots <- with(mam_subset, table(plotID, taxonID))
M_siteplotnights <- with(mam_subset, table(interaction(plotID, date), taxonID))

totalabund <- apply(M_sites,1,sum)
totalrich <- apply(M_sites>0,1,sum)
totalabundplot <- apply(M_siteplots,1,sum)
totalrichplot <- apply(M_siteplots>0,1,sum)

library(vegan)
library(ggplot2)

site_rare <- rarefy(x = M_sites, MARGIN=1, sample=min(totalabund), se=TRUE)
ggplot(data.frame(obs_richness=totalrich, rare_richness=site_rare[1,], se=site_rare[2,]),
       aes(x=obs_richness, y=rare_richness, ymin=rare_richness-se, ymax=rare_richness+se)) +
  geom_pointrange() + theme_minimal() +
  labs(x='Observed richness', y='Rarefied richness') +
  geom_abline(slope=1, intercept=0, linetype='dotted', color='red') +
  xlim(0,12) + ylim(0,12)
ggsave('figs/png/mammalrarefaction_bysite.png', height=5, width=5)

site_plot_rare <- rarefy(x = M_siteplots, MARGIN=1, sample=5, se=TRUE)


# # Calculate FDisp metric with and without intraspecific variation.
# 
# library(FD)
# mam_justtr <- mam_subset[,c('siteID','plotID','taxonID','hindfootLength','earLength','tailLength','totalLength','weight')]
# 
# # Means by species, ignoring sites
# # Convert to the proper format for the FD library.
# trmeans_sp <- ddply(mam_justtr, .(taxonID), colwise(mean, is.numeric, na.rm=TRUE))
# trmeans_sp[is.na(trmeans_sp)] <- NA
# trmeans_sp_mat <- as.matrix(trmeans_sp[,-1])
# dimnames(trmeans_sp_mat)[[1]] <- trmeans_sp[,1] 
# 
# # Get rid of all-NA rows
# useRows <- apply(trmeans_sp_mat, 1, function(x) !all(is.na(x)))
# useRows['GLVO'] <- FALSE # Including GLVO causes it to fail because it only has one trait recorded.
# useRows[c('GLVO','MIPE','PEMA')] <- FALSE # Exclude the same ones as from the other site.
# 
# gower_sp <- gowdis(trmeans_sp_mat[useRows,])
# fd_mam_bysite <- fdisp(d = gower_sp, a = M_sites[,useRows])
# 
# # Means by species accounting for AMONG site variation only
# trmeans_sitesp <- ddply(mam_justtr, .(siteID, taxonID), colwise(mean, is.numeric, na.rm=TRUE))
# trmeans_sitesp[is.na(trmeans_sitesp)] <- NA
# trmeans_sitesp_mat <- as.matrix(trmeans_sitesp[,-(1:2)])
# dimnames(trmeans_sitesp_mat)[[1]] <- paste(trmeans_sitesp[,1], trmeans_sitesp[,2], sep='.')
# 
# useRows <- apply(trmeans_sitesp_mat, 1, function(x) !all(is.na(x)))
# useRows[grep('GLVO|MIPE|PEMA', dimnames(trmeans_sitesp_mat)[[1]])] <- FALSE # Including GLVO causes it to fail because it only has one trait recorded.
# 
# gower_sitesp <- gowdis(trmeans_sitesp_mat[useRows, ])
# 
# # Make a new community matrix.
# M_sitesp <- with(mam_subset, table(siteID, interaction(siteID,taxonID)))
# M_sitesp <- M_sitesp[, dimnames(trmeans_sitesp_mat)[[1]]] # Match order to the trait matrix.
# 
# fd_mam_bysitesp <- fdisp(d = gower_sitesp, a=M_sitesp[, useRows])


# Now include ALL individuals, ignoring both species and site.
# Might be a problem due to the high number of NAs


########################################################
# Do all these FD metrics on the same individuals: ones that have at least weight AND hind foot length recorded.

# Do the subsetting by unique individual
# Subsetting
#mam_capture_subset <- subset(mam_capture, !siteID %in% c('DSNY')) not needed any more 21 July
mam_capture_subset <- mam_capture

# Take the median value of the traits for each individual if it was recaptured
library(plyr)
mam_noID <- mam_capture_subset[mam_capture_subset$individualID=='', c('siteID','taxonID','individualID','hindfootLength','earLength','tailLength','totalLength','weight')]
mam_byindiv <- ddply(subset(mam_capture_subset, individualID != ''),
                     .(siteID, taxonID, individualID), 
                     colwise(median, .(hindfootLength, earLength, tailLength, totalLength, weight), na.rm=TRUE))
mam_byindiv <-rbind(mam_byindiv, mam_noID)



library(FD)

mam_subset2 <- subset(mam_byindiv, !is.na(weight) & !is.na(hindfootLength)) # Retain 3989 rows
mam_justtr <- mam_subset2[,c('siteID','taxonID','hindfootLength','earLength','tailLength','totalLength','weight')]

# Using single mean for each species across all sites
trmeans_sp <- ddply(mam_justtr, .(taxonID), colwise(mean, is.numeric, na.rm=TRUE))
trmeans_sp[is.na(trmeans_sp)] <- NA
trmeans_sp_mat <- as.matrix(trmeans_sp[,-1])
dimnames(trmeans_sp_mat)[[1]] <- trmeans_sp[,1] 

# Using different mean for each species at each site
trmeans_sitesp <- ddply(mam_justtr, .(siteID, taxonID), colwise(mean, is.numeric, na.rm=TRUE))
trmeans_sitesp[is.na(trmeans_sitesp)] <- NA
trmeans_sitesp_mat <- as.matrix(trmeans_sitesp[,-(1:2)])
dimnames(trmeans_sitesp_mat)[[1]] <- paste(trmeans_sitesp[,1], trmeans_sitesp[,2], sep='.')

# Using all individuals' values
trvals_indiv <- as.matrix(mam_justtr[,c('hindfootLength','earLength','tailLength','totalLength','weight')])
dimnames(trvals_indiv)[[1]] <- paste(mam_justtr$siteID, mam_justtr$taxonID, 1:nrow(mam_justtr), sep='.')

# Calculate gower distance matrices (must not contain NAs)
gower_sp <- gowdis(trmeans_sp_mat)
gower_sitesp <- gowdis(trmeans_sitesp_mat)
gower_indiv <- gowdis(trvals_indiv)

# Create site by species matrices for each one and ensure labels all match and are in the correct order.
M_sites <- with(mam_subset2, table(siteID, taxonID))
M_sitesp <- with(mam_subset2, table(siteID, interaction(siteID,taxonID)))
M_sitesp <- M_sitesp[, dimnames(trmeans_sitesp_mat)[[1]]]
M_indiv <- with(mam_subset2, table(siteID, paste(mam_subset2$siteID,mam_subset2$taxonID,1:nrow(mam_subset2), sep='.')))
M_indiv <- M_indiv[, dimnames(trvals_indiv)[[1]]]

# Calculate FD.
fd_mam_bysite <- fdisp(d = gower_sp, a = M_sites)
fd_mam_bysitesp <- fdisp(d = gower_sitesp, a = M_sitesp)
fd_mam_byindiv <- fdisp(d = gower_indiv, a = M_indiv)


# Make a plot without the "WOOD" site that has 0 diversity

fd_mam_df <- data.frame(fd_spmeans = fd_mam_bysite$FDis, 
                        fd_spsitemeans = fd_mam_bysitesp$FDis,
                        fd_indivvalues = fd_mam_byindiv$FDis)

library(ggplot2)
library(GGally)

points_abline <- function(data, mapping, ...){
  p <- ggplot(data = data, mapping = mapping) + 
    geom_point() + 
    geom_abline(slope=1, intercept=0, linetype='dotted', color='red')
  p
}

png('figs/png/mammal_FD_methodscomparison.png', height=8, width=8, res=300, units = 'in')

ggpairs(fd_mam_df, 
        lower=list(continuous=points_abline),
        diag=list(continuous=wrap('barDiag', bins=8)),
        columnLabels = c('Global species\nmeans', 'Site-specific\nspecies means', 'Individual trait\nvalues')) +
  theme_bw()

dev.off()