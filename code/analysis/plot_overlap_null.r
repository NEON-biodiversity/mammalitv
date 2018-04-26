# Plot overlap stats from null models
# Author: QDR
# Project: NEON ITV
# Created: 30 Sep 2016
# Last modified: 20 Apr 2018

# Modified 20 Apr 2018: Specify paths to data on HPCC.
data_path <- '/mnt/research/neon'

# Load data
library(cati)
library(dplyr)
library(ggplot2)
library(lubridate)
library(reshape2)
source('code/analysis/T_statistics/tstat2longform.r')
source('code/analysis/densityoverlap.r')
source('~/GitHub/NEON/code/bioclimnames.r')

load(file.path(data_path, 'MS1_RodentOverlap/R_data/overlapstat/overlap_allyears_bysite.r'))
load(file.path(data_path, 'MS1_RodentOverlap/R_data/overlapstat/overlap_2015_bysite.r'))

load(file.path(data_path, 'final_data/allorganismal_latest.r'))
load(file.path(data_path, 'final_data/allsiteplot_latest.r'))
load(file.path(data_path, 'MS1_RodentOverlap/R_data/mammalPDbyplotobject.r'))

mammalTax <- read.csv(file.path(data_path, 'raw_data/NEON_mam_taxonomy.csv'))
mammalTraits <- read.csv(file.path(data_path, 'external_data/final_external_data/NEON_miscmammaltraits.csv'))

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

spatialstat <- read.csv(file.path(data_path, 'external_data/final_external_data/NEON_spatial_stats_highres.csv'))
spstatmeans <- neonplotdata %>% cbind(spatialstat) %>% group_by(siteID) %>% summarize(ruggedness=mean(tri, na.rm=TRUE))

chao <- function(x) {
  xcomm <- table(x$taxonID)
  S_obs <- length(xcomm)
  f1 <- sum(xcomm == 1)
  f2 <- sum(xcomm == 2)
  return(data.frame(chao1 = S_obs + (f1 * (f1 - 1)) / (2 * (f2 + 1))))
}

chao1site <- mam_capture %>% group_by(siteID) %>% do(chao(.))
chao15site <- mam_capture %>% filter(year==2015) %>% group_by(siteID) %>% do(chao(.))
mammalPDsite <- mammalPD %>% mutate(siteID = substr(plotID,1,4)) %>% group_by(siteID) %>% summarize(mpd_z = mean(mpd.obs.z2015, na.rm=T), mntd_z = mean(mntd.obs.z2015, na.rm=T))


ovlses_site_all <- overlap2longform_ses(overlap_allyears_bysite) %>% left_join(neonsitedata) %>% left_join(spstatmeans) %>% left_join(mammalPDsite) %>% left_join(chao1site) %>% filter(chao1 > 1)

ovlses_site_2015 <- overlap2longform_ses(overlap_2015_bysite) %>% left_join(neonsitedata) %>% left_join(spstatmeans) %>% left_join(mammalPDsite) %>% left_join(chao15site) %>% filter(chao1 > 1, trait == 'logweight', !siteID %in% c('DELA','HEAL','DSNY'))

ggplot(ovlses_site_all %>% filter(trait=='logweight', !siteID %in% c('DELA','HEAL','DSNY')), aes(x=bio1, y=overlap)) + geom_point()
ggplot(ovlses_site_all %>% filter(trait=='logweight', !siteID %in% c('DELA','HEAL','DSNY')), aes(x=bio1, y=overlap_ses)) + geom_point()
ggplot(ovlses_site_all %>% filter(trait=='logweight', !siteID %in% c('DELA','HEAL','DSNY')), aes(x=chao1, y=overlap)) + geom_point()
ggplot(ovlses_site_all %>% filter(trait=='logweight', !siteID %in% c('DELA','HEAL','DSNY')), aes(x=chao1, y=overlap_ses)) + geom_point()


# Scatterplots with null dists. -------------------------------------------

posestemp <- ggplot(ovlses_site_2015, aes(x=bio1)) + 
  geom_segment(aes(y = ci_min, yend = ci_max, xend=bio1), alpha = 0.5, size = 1.5, color = 'skyblue') +
  geom_point(aes(y = overlap_ses), size = 1.5) +
  labs(y = 'Overlap Statistic SES', x = parse(text=bioclimnames[1])) +
  theme_minimal() + theme(panel.border = element_rect(fill='transparent')) +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'midnightblue', size=0.9) +
  qSubtitle('Overlap statistics for 2015 NEON mammals versus MAT', 'Standardized Effect Sizes')

potemp <- ggplot(ovlses_site_2015, aes(x=bio1)) + 
  geom_segment(aes(y = ci_min_raw, yend = ci_max_raw, xend=bio1), alpha = 0.5, size = 1.5, color = 'skyblue') +
  geom_point(aes(y = overlap), size = 1.5) +
  labs(y = 'Overlap Statistic', x = parse(text=bioclimnames[1])) +
  theme_minimal() + theme(panel.border = element_rect(fill='transparent')) +
  #geom_hline(yintercept = 0, linetype = 'dotted', color = 'midnightblue', size=0.9) +
  qSubtitle('Overlap statistics for 2015 NEON mammals versus MAT', 'Raw Values')

poseschao <- ggplot(ovlses_site_2015, aes(x=chao1)) + 
  geom_segment(aes(y = ci_min, yend = ci_max, xend=chao1), alpha = 0.5, size = 1.5, color = 'skyblue') +
  geom_point(aes(y = overlap_ses), size = 1.5) +
  labs(y = 'Overlap Statistic SES', x = 'Chao1 Species Richness') +
  theme_minimal() + theme(panel.border = element_rect(fill='transparent')) +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'midnightblue', size=0.9) +
  qSubtitle('Overlap statistics for 2015 NEON mammals versus richness', 'Standardized Effect Sizes')

pochao <- ggplot(ovlses_site_2015, aes(x=chao1)) + 
  geom_segment(aes(y = ci_min_raw, yend = ci_max_raw, xend=chao1), alpha = 0.5, size = 1.5, color = 'skyblue') +
  geom_point(aes(y = overlap), size = 1.5) +
  labs(y = 'Overlap Statistic', x = 'Chao1 Species Richness') +
  theme_minimal() + theme(panel.border = element_rect(fill='transparent')) +
  #geom_hline(yintercept = 0, linetype = 'dotted', color = 'midnightblue', size=0.9) +
  qSubtitle('Overlap statistics for 2015 NEON mammals versus richness', 'Raw Values')

posesmpd <- ggplot(ovlses_site_2015, aes(x=mpd_z)) + 
  geom_segment(aes(y = ci_min, yend = ci_max, xend=mpd_z), alpha = 0.5, size = 1.5, color = 'skyblue') +
  geom_point(aes(y = overlap_ses), size = 1.5) +
  labs(y = 'Overlap Statistic SES', x = 'Pairwise Phy Dist SES') +
  theme_minimal() + theme(panel.border = element_rect(fill='transparent')) +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'midnightblue', size=0.9) +
  qSubtitle('Overlap statistics for 2015 NEON mammals versus MPD', 'Standardized Effect Sizes')

pompd <- ggplot(ovlses_site_2015, aes(x=mpd_z)) + 
  geom_segment(aes(y = ci_min_raw, yend = ci_max_raw, xend=mpd_z), alpha = 0.5, size = 1.5, color = 'skyblue') +
  geom_point(aes(y = overlap), size = 1.5) +
  labs(y = 'Overlap Statistic', x = 'Pairwise Phy Dist SES') +
  theme_minimal() + theme(panel.border = element_rect(fill='transparent')) +
  #geom_hline(yintercept = 0, linetype = 'dotted', color = 'midnightblue', size=0.9) +
  qSubtitle('Overlap statistics for 2015 NEON mammals versus MPD', 'Raw Values')