# Load new O-stat data (using no opportunistic species, using harmonic means, and using the two different types of null model--original and Mason's)
# QDR 12 Jan 2016

# Modified 26 Apr 2018: Specify path to HPCC at top
# "Fork" 04 May: yet another new null model, swapping means but keeping individuals' deviations. See whether this is interesting.
# Modified 16 Jan: Add OAES in (the 2015 site that was inadvertently left out of the older data objects)

data_path <- '/mnt/research/neon'

library(dplyr)
library(ggplot2)
library(lubridate)
library(reshape2)
library(extrafont)

fp <- file.path(data_path, 'MS1_RodentOverlap/R_data/overlapstat/12Jan')
source('code/analysis/ostat2longform.r')
source('~/GitHub/NEON/code/bioclimnames.r')

# New null
load(file.path(fp, "Ostats_bysite2015_harmoniciucn_shufflewt.r"))
Ostats_bysite2015iucn_newnull <- Ostats_bysite2015iucn

# New"er" null
load(file.path(data_path, 'MS1_RodentOverlap/R_data/overlapstat/03May/Ostats_bysite2015_harmoniciucn_swapmeans.r'))
Ostats_bysite2015_iucn_swapnull <- Ostats_bysite2015iucn

# Original null
load(file.path(fp, "Ostats_bysite2015_harmoniciucn_orignull.r"))

o2015 <- ostat2longform(Ostats_bysite2015iucn)
o2015n <- ostat2longform(Ostats_bysite2015iucn_newnulls)
o2015sn <- ostat2longform(Ostats_bysite2015_iucn_swapnull)

iucn <- read.csv(file.path(data_path, 'external_data/final_external_data/IUCN_mammal_ranges.csv'))
load(file.path(data_path, 'final_data/allorganismal_latest.r'))
load(file.path(data_path, 'final_data/allsiteplot_latest.r'))

mammalTax <- read.csv(file.path(data_path, 'final_data/mammals/NEON_mam_taxonomy.csv'))
mammalTraits <- read.csv(file.path(data_path, 'external_data/final_external_data/NEON_miscmammaltraits.csv'))

# We need new richness and phylogenetic diversity calculations, because the old ones are based on the wrong species list. Do everything from scratch here?

# 12 Jan: calculated new mammal phylogenetic diversity. Updated plot object.
load(file.path(data_path, 'MS1_RodentOverlap/R_data/mammalPDobject.r'))

spatialstat <- read.csv(file.path(data_path, 'external_data/final_external_data/NEON_spatial_stats_highres.csv'))
spstatmeans <- neonplotdata %>% cbind(spatialstat) %>% group_by(siteID) %>% summarize(ruggedness=mean(tri, na.rm=TRUE))

heterodf <- read.csv(file.path(data_path, 'MS1_RodentOverlap/R_data/heterogeneity.csv'), stringsAsFactors = FALSE)
mammalPDsite <- mammalPD %>% rename(mpd_z = mpd.obs.z2015, mntd_z = mntd.obs.z2015) %>% dplyr::select(siteID, mpd_z, mntd_z)

# Recalc richness estimators
library(iNEXT)
mam_forrich <- mutate(mam_capture, 
                      individualandtag = pmin(as.character(individualID), as.character(tagID), na.rm=TRUE),
                      year = year(date)) %>%
  filter(order == 'Rodentia', taxonProtocolCategory == 'target', year < 2016)

# Get rid of poorly identified individuals (those marked sp.) because they should not count for our sampling . . . they are spurious singletons.

mammaltables <- mam_forrich %>% filter(!siteID %in% c('DSNY','DELA','HEAL'), !grepl('sp\\.', mam_forrich$scientificName)) %>% group_by(siteID) %>% do(t = table(.$taxonID))

# Make list of abundances into a species X site matrix (species must be rows, sites must be columns)
mammaltaxalist <- unique(unlist(lapply(mammaltables$t, names)))
mammalmat <- sapply(mammaltables$t, function(x) x[match(mammaltaxalist, names(x))])
mammalmat[is.na(mammalmat)] <- 0
dimnames(mammalmat) <- list(mammaltaxalist, mammaltables$siteID)

# With list
mamx <- lapply(mammaltables$t, as.numeric)
names(mamx) <- mammaltables$siteID

set.seed(46545)
default <- iNEXT(x=mamx, q=0, datatype='abundance', size = c(5,10,50,100,500,1000,2000,3000,4000))

# Asymptotic estimator compared to Chao1

chao <- function(x) {
  xcomm <- table(x$taxonID)
  S_obs <- length(xcomm)
  f1 <- sum(xcomm == 1)
  f2 <- sum(xcomm == 2)
  return(data.frame(chao1 = S_obs + (f1 * (f1 - 1)) / (2 * (f2 + 1))))
}



asyrich <- default$AsyEst %>% filter(Diversity == 'Species richness') 
asyrich$Site <- factor(asyrich$Site, levels=asyrich$Site[order(asyrich$Observed)])

chao1site <- mam_forrich %>% filter(!siteID %in% c('DSNY','DELA','HEAL'), !grepl('sp\\.', mam_forrich$scientificName)) %>% group_by(siteID) %>% do(chao(.))
chao1site$siteID <-factor(asyrich$Site, levels=asyrich$Site[order(asyrich$Observed)])

richnessdat <- full_join(asyrich %>% dplyr::select(Site, Observed, Estimator) %>% rename(siteID=Site), chao1site)

# Create new mam_capture_sitemerge.

mammalGuilds <- mammalTraits %>% 
  mutate(scientificName = paste(Genus, Species)) %>%
  dplyr::select(scientificName, Pineda_Main_food)

mammalGuilds <- rbind(mammalGuilds, data.frame(scientificName=c("Dipodomys sp.", "Glaucomys sp.", "Microtus sp.", "Neotoma sp.", 
                                                                "Perognathus sp.", "Peromyscus sp.", "Reithrodontomys sp."),
                                               Pineda_Main_food=c('Generalist','Generalist','Herbivore','Herbivore','Generalist','Granivore','Generalist')))


mam_forrich <- left_join(mam_forrich, mammalGuilds[-c(89,118,62,63),], by = 'scientificName')

mam_noID <- mam_forrich[mam_forrich$individualandtag=='', c('year', 'siteID','plotID', 'taxonID','family', 'individualandtag','hindfootLength','earLength','tailLength','totalLength','weight','sex','lifeStage','Pineda_Main_food')]

mam_grp <- filter(mam_forrich, individualandtag != '') %>%
  group_by(year, siteID, plotID, taxonID, family, individualandtag)

mam_byindiv <- summarize_at(mam_grp, vars(hindfootLength, earLength, tailLength, totalLength, weight), median, na.rm=TRUE)
mam_byindiv_class <- do(mam_grp, data.frame(sex = .$sex[1], lifeStage = names(sort(table(.$lifeStage),decreasing=TRUE))[1], Pineda_Main_food=.$Pineda_Main_food[1]))

mam_byindiv <- cbind(as.data.frame(mam_byindiv), as.data.frame(mam_byindiv_class[,c('sex', 'lifeStage','Pineda_Main_food')]))                  
mam_byindiv <-rbind(mam_byindiv, mam_noID)

mam_capture_sitemerge <- left_join(mam_byindiv, filter(neonplotdata, subtype == 'mammalGrid'))


# Join err'thang togetha

o2015 <- o2015 %>% rename(siteID=site) %>% left_join(neonsitedata) %>% left_join(spstatmeans) %>% left_join(mammalPDsite) %>% left_join(richnessdat) %>% left_join(heterodf) 
o2015n <- o2015n %>% rename(siteID=site) %>% left_join(neonsitedata) %>% left_join(spstatmeans) %>% left_join(mammalPDsite) %>% left_join(richnessdat) %>% left_join(heterodf)
o2015sn <- o2015sn %>% rename(siteID=site) %>% left_join(neonsitedata) %>% left_join(spstatmeans) %>% left_join(mammalPDsite) %>% left_join(richnessdat) %>% left_join(heterodf)


# Set plot parameters
theme_john <- theme_bw() + theme(panel.grid = element_blank(), 
                                 axis.text = element_text(size = 12), 
                                 axis.title = element_text(size = 18),
                                 text = element_text(family = 'Helvetica'))

# Add some additional parameters to get rid of the axis ticks and numbers if you want.
theme_noaxisnumbers <- theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks.x = element_blank(), axis.ticks.y = element_blank())


# Regressions -------------------------------------------------------------

library(betareg)
library(MuMIn)

reglocal <- betareg(ostat_norm ~ bio1 + chao1 + mpd_z + pc1_productivityheterogeneity + pc2_topographyheterogeneity, link = 'logit', data = o2015 %>% filter(trait == 'logweight'), na.action = 'na.pass')

summary(reglocal)
confint(reglocal)
dredge(reglocal)

reglocalbest <- betareg(ostat_norm ~ bio1 + chao1, link = 'logit', data = o2015 %>% filter(trait == 'logweight'), na.action = 'na.pass')

summary(reglocalbest)
confint(reglocalbest)

# check for multicollinearity

library(usdm)
pred_df <- o2015 %>% filter(trait == 'logweight') %>% dplyr::select(bio1, chao1, mpd_z, pc1_productivityheterogeneity, pc2_topographyheterogeneity)
pairs(pred_df)
cor(pred_df)
vif(pred_df) # Looks OK.

# Look at distributions of (local) effect sizes for both null models.

o2015$ostat_norm_localnull_ses[o2015$trait == 'logweight'] # Original null.
o2015n$ostat_norm_localnull_ses[o2015n$trait == 'logweight'] # "Mason" null.
o2015sn$ostat_norm_localnull_ses[o2015sn$trait == 'logweight'] # "Swap means" null.

o2015sn[o2015sn$trait == 'logweight', c('siteID','ostat_norm_localnull_ses')] # "Swap means" null.


quantile(o2015$ostat_norm_localnull_ses[o2015$trait == 'logweight'], probs = c(0.5, 0.025, 0.975)) # Median and 95% envelope
quantile(o2015sn$ostat_norm_localnull_ses[o2015sn$trait == 'logweight'], probs = c(0.5, 0.025, 0.975), na.rm=T) # Median and 95% envelope


o2015 <- o2015 %>%
  mutate(local_significant = ostat_norm_localnull_ses < ostat_norm_localnull_seslower | ostat_norm_localnull_ses > ostat_norm_localnull_sesupper,
         reg_significant = ostat_norm_regnull_ses < ostat_norm_regnull_seslower | ostat_norm_regnull_ses > ostat_norm_regnull_sesupper)
o2015n <- o2015n %>%
  mutate(local_significant = ostat_norm_localnull_ses < ostat_norm_localnull_seslower | ostat_norm_localnull_ses > ostat_norm_localnull_sesupper,
         reg_significant = ostat_norm_regnull_ses < ostat_norm_regnull_seslower | ostat_norm_regnull_ses > ostat_norm_regnull_sesupper)

table(o2015$local_significant[o2015$trait == 'logweight']) # 5 of 23 in neutral zone, 18 of 23 are partitioned.
table(o2015n$local_significant[o2015n$trait == 'logweight']) # All but 1 are in neutral zone. No evidence that commonness has anything to do with it, as it does for plants, which Mason showed.

# New plots ---------------------------------------------------------------

reglocalbio <- betareg(ostat_norm ~ bio1, link = 'logit', data = o2015 %>% filter(trait == 'logweight'))
reglocalchao <- betareg(ostat_norm ~ chao1, link = 'logit', data = o2015 %>% filter(trait == 'logweight'))

fx <- function(x, b0, b1) 1/(1 + exp(-(b0 + b1 * x)))
fx_loglog <- function(x, b0, b1) exp(-exp(-(b0 + b1 * x)))
tempco <- summary(reglocalbio)$coeff$mean[,1]
chaoco <- summary(reglocalchao)$coeff$mean[,1]
csc <- scale_color_manual(values = c('gray75','black'))
tempstring <- "paste(\"Mean Annual Temperature (\",degree,\"C)\")" 


