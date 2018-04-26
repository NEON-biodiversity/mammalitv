# New mammal t-stat plots. See below for details.
# Author: QDR
# Project: NEON ITV
# Created: 27 Sep 2016
# Last Modified: 28 Sep 2016

# Q tweaked the T-stat calculation (see tstat_weightedcalc.r) to weight by relative abundance where appropriate.
# also, all guilds are split out individually. We did not include insectivorous rodents or shrews.

# Modified 28 Sep: Also look at raw data (not just SES)

# Load and clean data -----------------------------------------------------


library(cati)
library(dplyr)
library(ggplot2)
library(lubridate)
library(reshape2)
source('~/GitHub/NEON/code/analysis/tstat2longform.r')
source('~/qutil.r')
source('~/GitHub/NEON/code/bioclimnames.r')

load('~/GitHub/NEON/allorganismal_latest.r')
load('~/GitHub/NEON/allsiteplot_latest.r')
load('C:/Users/Q/Dropbox/neon/code/mammalPDbyplotobject.r')

load('C:/Users/Q/Dropbox/neon/code/tstatweighted/tstats_2015_gen_iucn.r')
load('C:/Users/Q/Dropbox/neon/code/tstatweighted/tstats_2015_gg_iucn.r')
load('C:/Users/Q/Dropbox/neon/code/tstatweighted/tstats_2015_ggh_iucn.r')
load('C:/Users/Q/Dropbox/neon/code/tstatweighted/tstats_2015_gh_iucn.r')
load('C:/Users/Q/Dropbox/neon/code/tstatweighted/tstats_2015_gran_iucn.r')
load('C:/Users/Q/Dropbox/neon/code/tstatweighted/tstats_2015_herb_iucn.r')

mammalTax <- read.csv('C:/Users/Q/Dropbox/neon/protocols/taxonomy/NEON_mam_taxonomy.csv')
mammalTraits <- read.csv('C:/Users/Q/Dropbox/neon/data/external_datasets/NEON_miscmammaltraits.csv')


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

spatialstat <- read.csv('C:/Users/Q/Dropbox/neon/data/external_datasets/NEON_spatial_stats_highres.csv')
neonplotdata <- cbind(neonplotdata, spatialstat) %>% rename(ruggedness = tri)

tmgg <- tstat2longform_ses(tstats_2015_gg_iucn, bysite = FALSE) %>% 
  left_join(filter(neonplotdata, subtype == 'mammalGrid'))
tmgh <- tstat2longform_ses(tstats_2015_gh_iucn, bysite = FALSE) %>% 
  left_join(filter(neonplotdata, subtype == 'mammalGrid'))
tmggh <- tstat2longform_ses(tstats_2015_ggh_iucn, bysite = FALSE) %>% 
  left_join(filter(neonplotdata, subtype == 'mammalGrid'))
tmgen <- tstat2longform_ses(tstats_2015_gen_iucn, bysite = FALSE) %>% 
  left_join(filter(neonplotdata, subtype == 'mammalGrid'))
tmgran <- tstat2longform_ses(tstats_2015_gran_iucn, bysite = FALSE) %>% 
  left_join(filter(neonplotdata, subtype == 'mammalGrid'))
tmherb <- tstat2longform_ses(tstats_2015_herb_iucn, bysite = FALSE) %>% 
  left_join(filter(neonplotdata, subtype == 'mammalGrid'))

# Modification: use raw data not effect sizes
tmgg <- tstat2longform(tstats_2015_gg_iucn, bysite = FALSE) %>% 
  left_join(filter(neonplotdata, subtype == 'mammalGrid'))
tmgh <- tstat2longform(tstats_2015_gh_iucn, bysite = FALSE) %>% 
  left_join(filter(neonplotdata, subtype == 'mammalGrid'))
tmggh <- tstat2longform(tstats_2015_ggh_iucn, bysite = FALSE) %>% 
  left_join(filter(neonplotdata, subtype == 'mammalGrid'))
tmgen <- tstat2longform(tstats_2015_gen_iucn, bysite = FALSE) %>% 
  left_join(filter(neonplotdata, subtype == 'mammalGrid'))
tmgran <- tstat2longform(tstats_2015_gran_iucn, bysite = FALSE) %>% 
  left_join(filter(neonplotdata, subtype == 'mammalGrid'))
tmherb <- tstat2longform(tstats_2015_herb_iucn, bysite = FALSE) %>% 
  left_join(filter(neonplotdata, subtype == 'mammalGrid'))

# Mammals have abundance so we use the Chao1 estimator.

chao <- function(x) {
  xcomm <- table(x$taxonID)
  S_obs <- length(xcomm)
  f1 <- sum(xcomm == 1)
  f2 <- sum(xcomm == 2)
  return(data.frame(chao1 = S_obs + (f1 * (f1 - 1)) / (2 * (f2 + 1))))
}

chao1gg15 <- mam_capture %>% filter(year(date) == 2015, Pineda_Main_food %in% c('Generalist','Granivore')) %>%
  group_by(plotID) %>% do(chao(.))
chao1gh15 <- mam_capture %>% filter(year(date) == 2015, Pineda_Main_food %in% c('Generalist','Herbivore')) %>%
  group_by(plotID) %>% do(chao(.))
chao1ggh15 <- mam_capture %>% filter(year(date) == 2015, Pineda_Main_food %in% c('Generalist','Granivore','Herbivore')) %>%
  group_by(plotID) %>% do(chao(.))
chao1gen15 <- mam_capture %>% filter(year(date) == 2015, Pineda_Main_food %in% c('Generalist')) %>%
  group_by(plotID) %>% do(chao(.))
chao1gran15 <- mam_capture %>% filter(year(date) == 2015, Pineda_Main_food %in% c('Granivore')) %>%
  group_by(plotID) %>% do(chao(.))
chao1herb15 <- mam_capture %>% filter(year(date) == 2015, Pineda_Main_food %in% c('Herbivore')) %>%
  group_by(plotID) %>% do(chao(.))

# Combine data frames, getting rid of the one species communities
tmgg <- tmgg %>% left_join(chao1gg15) %>% left_join(mammalPD) %>% filter(chao1 > 1)
tmgh <- tmgh %>% left_join(chao1gh15) %>% left_join(mammalPD) %>% filter(chao1 > 1)
tmggh <- tmggh %>% left_join(chao1ggh15) %>% left_join(mammalPD) %>% filter(chao1 > 1)
tmgen <- tmgen %>% left_join(chao1gen15) %>% left_join(mammalPD) %>% filter(chao1 > 1)
tmgran <- tmgran %>% left_join(chao1gran15) %>% left_join(mammalPD) %>% filter(chao1 > 1)
tmherb <- tmherb %>% left_join(chao1herb15) %>% left_join(mammalPD) %>% filter(chao1 > 1)



# Scatterplot -------------------------------------------------------------


pttempgg <- ggplot(subset(tmgg, trait == 'logweight'), aes(x=bio1)) + 
  facet_wrap( ~ stat, scales = 'free') +
  #geom_segment(aes(y = ci_min, yend = ci_max, xend=bio1), alpha = 0.5, size = 1.5, color = 'skyblue') +
  geom_point(aes(y = t), size = 1.5) +
  labs(y = 'T-statistic', x = parse(text=bioclimnames[1])) +
  theme_minimal() + theme(panel.border = element_rect(fill='transparent')) +
  #geom_hline(yintercept = 0, linetype = 'dotted', color = 'midnightblue', size=0.9) +
  qSubtitle('T-statistics for 2015 NEON mammals versus MAT', 'Generalists + Granivores')

ptmpdgg <- ggplot(subset(tmgg, trait == 'logweight'), aes(x=mpd.obs.z2015)) + 
  facet_grid(~ stat, scales = 'free') +
  #geom_segment(aes(y = ci_min, yend = ci_max, xend=mpd.obs.z2015), alpha = 0.5, size = 1.5, color = 'skyblue') +
  geom_point(aes(y = t), size = 1.5) +
  labs(y = 'T-statistic', x = 'SES of mean pairwise distance') +
  theme_minimal() + theme(panel.border = element_rect(fill='transparent')) +
  #geom_hline(yintercept = 0, linetype = 'dotted', color = 'midnightblue', size=0.9) +
  qSubtitle('T-statistics for 2015 NEON mammals versus MPD', 'Generalists + Granivores')

ptchaogg <- ggplot(subset(tmgg, trait == 'logweight'), aes(x=chao1)) + 
  facet_grid(~ stat, scales = 'free') +
  #geom_segment(aes(y = ci_min, yend = ci_max, xend=chao1), alpha = 0.5, size = 1.5, color = 'skyblue') +
  geom_point(aes(y = t), size = 1.5) +
  labs(y = 'T-statistic', x = 'Richness (Chao1)') +
  theme_minimal() + theme(panel.border = element_rect(fill='transparent')) +
  #geom_hline(yintercept = 0, linetype = 'dotted', color = 'midnightblue', size=0.9) +
  qSubtitle('T-statistics for 2015 NEON mammals versus richness', 'Generalists + Granivores')

pttempgh <- ggplot(subset(tmgh, trait == 'logweight'), aes(x=bio1)) + 
  facet_wrap( ~ stat, scales = 'free') +
  #geom_segment(aes(y = ci_min, yend = ci_max, xend=bio1), alpha = 0.5, size = 1.5, color = 'skyblue') +
  geom_point(aes(y = t), size = 1.5) +
  labs(y = 'T-statistic', x = parse(text=bioclimnames[1])) +
  theme_minimal() + theme(panel.border = element_rect(fill='transparent')) +
  #geom_hline(yintercept = 0, linetype = 'dotted', color = 'midnightblue', size=0.9) +
  qSubtitle('T-statistics for 2015 NEON mammals versus MAT', 'Generalists + Herbivores')

ptmpdgh <- ggplot(subset(tmgh, trait == 'logweight'), aes(x=mpd.obs.z2015)) + 
  facet_grid(~ stat, scales = 'free') +
  #geom_segment(aes(y = ci_min, yend = ci_max, xend=mpd.obs.z2015), alpha = 0.5, size = 1.5, color = 'skyblue') +
  geom_point(aes(y = t), size = 1.5) +
  labs(y = 'T-statistic', x = 'SES of mean pairwise distance') +
  theme_minimal() + theme(panel.border = element_rect(fill='transparent')) +
  #geom_hline(yintercept = 0, linetype = 'dotted', color = 'midnightblue', size=0.9) +
  qSubtitle('T-statistics for 2015 NEON mammals versus MPD', 'Generalists + Herbivores')

ptchaogh <- ggplot(subset(tmgh, trait == 'logweight'), aes(x=chao1)) + 
  facet_grid(~ stat, scales = 'free') +
  #geom_segment(aes(y = ci_min, yend = ci_max, xend=chao1), alpha = 0.5, size = 1.5, color = 'skyblue') +
  geom_point(aes(y = t), size = 1.5) +
  labs(y = 'T-statistic', x = 'Richness (Chao1)') +
  theme_minimal() + theme(panel.border = element_rect(fill='transparent')) +
  #geom_hline(yintercept = 0, linetype = 'dotted', color = 'midnightblue', size=0.9) +
  qSubtitle('T-statistics for 2015 NEON mammals versus richness', 'Generalists + Herbivores')

pttempggh <- ggplot(subset(tmggh, trait == 'logweight'), aes(x=bio1)) + 
  facet_wrap( ~ stat, scales = 'free') +
  #geom_segment(aes(y = ci_min, yend = ci_max, xend=bio1), alpha = 0.5, size = 1.5, color = 'skyblue') +
  geom_point(aes(y = t), size = 1.5) +
  labs(y = 'T-statistic', x = parse(text=bioclimnames[1])) +
  theme_minimal() + theme(panel.border = element_rect(fill='transparent')) +
  #geom_hline(yintercept = 0, linetype = 'dotted', color = 'midnightblue', size=0.9) +
  qSubtitle('T-statistics for 2015 NEON mammals versus MAT', 'All Three Guilds')

ptmpdggh <- ggplot(subset(tmggh, trait == 'logweight'), aes(x=mpd.obs.z2015)) + 
  facet_grid(~ stat, scales = 'free') +
  #geom_segment(aes(y = ci_min, yend = ci_max, xend=mpd.obs.z2015), alpha = 0.5, size = 1.5, color = 'skyblue') +
  geom_point(aes(y = t), size = 1.5) +
  labs(y = 'T-statistic', x = 'SES of mean pairwise distance') +
  theme_minimal() + theme(panel.border = element_rect(fill='transparent')) +
  #geom_hline(yintercept = 0, linetype = 'dotted', color = 'midnightblue', size=0.9) +
  qSubtitle('T-statistics for 2015 NEON mammals versus MPD', 'All Three Guilds')

ptchaoggh <- ggplot(subset(tmggh, trait == 'logweight'), aes(x=chao1)) + 
  facet_grid(~ stat, scales = 'free') +
  #geom_segment(aes(y = ci_min, yend = ci_max, xend=chao1), alpha = 0.5, size = 1.5, color = 'skyblue') +
  geom_point(aes(y = t), size = 1.5) +
  labs(y = 'T-statistic', x = 'Richness (Chao1)') +
  theme_minimal() + theme(panel.border = element_rect(fill='transparent')) +
  #geom_hline(yintercept = 0, linetype = 'dotted', color = 'midnightblue', size=0.9) +
  qSubtitle('T-statistics for 2015 NEON mammals versus richness', 'All Three Guilds')

pttempgen <- ggplot(subset(tmgen, trait == 'logweight'), aes(x=bio1)) + 
  facet_wrap( ~ stat, scales = 'free') +
  #geom_segment(aes(y = ci_min, yend = ci_max, xend=bio1), alpha = 0.5, size = 1.5, color = 'skyblue') +
  geom_point(aes(y = t), size = 1.5) +
  labs(y = 'T-statistic', x = parse(text=bioclimnames[1])) +
  theme_minimal() + theme(panel.border = element_rect(fill='transparent')) +
  #geom_hline(yintercept = 0, linetype = 'dotted', color = 'midnightblue', size=0.9) +
  qSubtitle('T-statistics for 2015 NEON mammals versus MAT', 'Generalists only')

ptmpdgen <- ggplot(subset(tmgen, trait == 'logweight'), aes(x=mpd.obs.z2015)) + 
  facet_grid(~ stat, scales = 'free') +
  #geom_segment(aes(y = ci_min, yend = ci_max, xend=mpd.obs.z2015), alpha = 0.5, size = 1.5, color = 'skyblue') +
  geom_point(aes(y = t), size = 1.5) +
  labs(y = 'T-statistic', x = 'SES of mean pairwise distance') +
  theme_minimal() + theme(panel.border = element_rect(fill='transparent')) +
  #geom_hline(yintercept = 0, linetype = 'dotted', color = 'midnightblue', size=0.9) +
  qSubtitle('T-statistics for 2015 NEON mammals versus MPD', 'Generalists only')

ptchaogen <- ggplot(subset(tmgen, trait == 'logweight'), aes(x=chao1)) + 
  facet_grid(~ stat, scales = 'free') +
  #geom_segment(aes(y = ci_min, yend = ci_max, xend=chao1), alpha = 0.5, size = 1.5, color = 'skyblue') +
  geom_point(aes(y = t), size = 1.5) +
  labs(y = 'T-statistic', x = 'Richness (Chao1)') +
  theme_minimal() + theme(panel.border = element_rect(fill='transparent')) +
  #geom_hline(yintercept = 0, linetype = 'dotted', color = 'midnightblue', size=0.9) +
  qSubtitle('T-statistics for 2015 NEON mammals versus richness', 'Generalists only')

pttempgran <- ggplot(subset(tmgran, trait == 'logweight'), aes(x=bio1)) + 
  facet_wrap( ~ stat, scales = 'free') +
  #geom_segment(aes(y = ci_min, yend = ci_max, xend=bio1), alpha = 0.5, size = 1.5, color = 'skyblue') +
  geom_point(aes(y = t), size = 1.5) +
  labs(y = 'T-statistic', x = parse(text=bioclimnames[1])) +
  theme_minimal() + theme(panel.border = element_rect(fill='transparent')) +
  #geom_hline(yintercept = 0, linetype = 'dotted', color = 'midnightblue', size=0.9) +
  qSubtitle('T-statistics for 2015 NEON mammals versus MAT', 'Granivores only')

ptmpdgran <- ggplot(subset(tmgran, trait == 'logweight'), aes(x=mpd.obs.z2015)) + 
  facet_grid(~ stat, scales = 'free') +
  #geom_segment(aes(y = ci_min, yend = ci_max, xend=mpd.obs.z2015), alpha = 0.5, size = 1.5, color = 'skyblue') +
  geom_point(aes(y = t), size = 1.5) +
  labs(y = 'T-statistic', x = 'SES of mean pairwise distance') +
  theme_minimal() + theme(panel.border = element_rect(fill='transparent')) +
  #geom_hline(yintercept = 0, linetype = 'dotted', color = 'midnightblue', size=0.9) +
  qSubtitle('T-statistics for 2015 NEON mammals versus MPD', 'Granivores only')

ptchaogran <- ggplot(subset(tmgran, trait == 'logweight'), aes(x=chao1)) + 
  facet_grid(~ stat, scales = 'free') +
  #geom_segment(aes(y = ci_min, yend = ci_max, xend=chao1), alpha = 0.5, size = 1.5, color = 'skyblue') +
  geom_point(aes(y = t), size = 1.5) +
  labs(y = 'T-statistic', x = 'Richness (Chao1)') +
  theme_minimal() + theme(panel.border = element_rect(fill='transparent')) +
  #geom_hline(yintercept = 0, linetype = 'dotted', color = 'midnightblue', size=0.9) +
  qSubtitle('T-statistics for 2015 NEON mammals versus richness', 'Granivores only')

pttempherb <- ggplot(subset(tmherb, trait == 'logweight'), aes(x=bio1)) + 
  facet_wrap( ~ stat, scales = 'free') +
  #geom_segment(aes(y = ci_min, yend = ci_max, xend=bio1), alpha = 0.5, size = 1.5, color = 'skyblue') +
  geom_point(aes(y = t), size = 1.5) +
  labs(y = 'T-statistic', x = parse(text=bioclimnames[1])) +
  theme_minimal() + theme(panel.border = element_rect(fill='transparent')) +
  #geom_hline(yintercept = 0, linetype = 'dotted', color = 'midnightblue', size=0.9) +
  qSubtitle('T-statistics for 2015 NEON mammals versus MAT', 'Herbivores only')

ptmpdherb <- ggplot(subset(tmherb, trait == 'logweight'), aes(x=mpd.obs.z2015)) + 
  facet_grid(~ stat, scales = 'free') +
  #geom_segment(aes(y = ci_min, yend = ci_max, xend=mpd.obs.z2015), alpha = 0.5, size = 1.5, color = 'skyblue') +
  geom_point(aes(y = t), size = 1.5) +
  labs(y = 'T-statistic', x = 'SES of mean pairwise distance') +
  theme_minimal() + theme(panel.border = element_rect(fill='transparent')) +
  #geom_hline(yintercept = 0, linetype = 'dotted', color = 'midnightblue', size=0.9) +
  qSubtitle('T-statistics for 2015 NEON mammals versus MPD', 'Herbivores only')

ptchaoherb <- ggplot(subset(tmherb, trait == 'logweight'), aes(x=chao1)) + 
  facet_grid(~ stat, scales = 'free') +
  #geom_segment(aes(y = ci_min, yend = ci_max, xend=chao1), alpha = 0.5, size = 1.5, color = 'skyblue') +
  geom_point(aes(y = t), size = 1.5) +
  labs(y = 'T-statistic', x = 'Richness (Chao1)') +
  theme_minimal() + theme(panel.border = element_rect(fill='transparent')) +
  #geom_hline(yintercept = 0, linetype = 'dotted', color = 'midnightblue', size=0.9) +
  qSubtitle('T-statistics for 2015 NEON mammals versus richness', 'Herbivores only')


pdf('figs/pdf/tstat_guild_scatterplots_rawvalues.pdf', height = 6, width = 12)
pttempgg
ptmpdgg
ptchaogg

pttempgh
ptmpdgh
ptchaogh

pttempggh
ptmpdggh
ptchaoggh

pttempgen
ptmpdgen
ptchaogen

pttempgran
ptmpdgran
ptchaogran

pttempherb
ptmpdherb
ptchaoherb
dev.off()

# Boxplot -----------------------------------------------------------------

bpgg <- ggplot(filter(tmgg, trait == 'logweight'), aes(x = stat, y = t)) +
  geom_jitter(height = 0, width = 0.25) + theme_bw() +
  ggtitle('Generalists + Granivores')
bpgh <- ggplot(filter(tmgh, trait == 'logweight'), aes(x = stat, y = t)) +
  geom_jitter(height = 0, width = 0.25) + theme_bw() +
  ggtitle('Generalists + Herbivores')
bpggh <- ggplot(filter(tmggh, trait == 'logweight'), aes(x = stat, y = t)) +
  geom_jitter(height = 0, width = 0.25) + theme_bw() +
  ggtitle('All three guilds')
bpgen <- ggplot(filter(tmgen, trait == 'logweight'), aes(x = stat, y = t)) +
  geom_jitter(height = 0, width = 0.25) + theme_bw() +
  ggtitle('Generalists only')
bpgran <- ggplot(filter(tmgran, trait == 'logweight'), aes(x = stat, y = t)) +
  geom_jitter(height = 0, width = 0.25) + theme_bw() +
  ggtitle('Granivores only')
bpherb <- ggplot(filter(tmherb, trait == 'logweight'), aes(x = stat, y = t)) +
  geom_jitter(height = 0, width = 0.25) + theme_bw() +
  ggtitle('Herbivores only')

pdf('figs/pdf/tstatboxplotsbyguild_rawvalues.pdf', height=6, width=6)
  bpgg
  bpgh
  bpggh
  bpgen
  bpgran
  bpherb
dev.off()


# Distribution plots ------------------------------------------------------


densgg <- ggplot(filter(tmgg, trait == 'logweight'), aes(color = stat, x = t)) +
  geom_density() + theme_bw() +
  ggtitle('Generalists + Granivores')
densgh <- ggplot(filter(tmgh, trait == 'logweight'), aes(color = stat, x = t)) +
  geom_density() + theme_bw() +
  ggtitle('Generalists + Herbivores')
densggh <- ggplot(filter(tmggh, trait == 'logweight'), aes(color = stat, x = t)) +
  geom_density() + theme_bw() +
  ggtitle('All three guilds')
densgen <- ggplot(filter(tmgen, trait == 'logweight'), aes(color = stat, x = t)) +
  geom_density() + theme_bw() +
  ggtitle('Generalists only')
densgran <- ggplot(filter(tmgran, trait == 'logweight'), aes(color = stat, x = t)) +
  geom_density() + theme_bw() +
  ggtitle('Granivores only')
densherb <- ggplot(filter(tmherb, trait == 'logweight'), aes(color = stat, x = t)) +
  geom_density() + theme_bw() +
  ggtitle('Herbivores only')

pdf('figs/pdf/tstatdensityplotsbyguild.pdf', height=6, width=6)
densgg
densgh
densggh
densgen
densgran
densherb
dev.off()



# Regression --------------------------------------------------------------

# Just do regular frequentist mixed model at the moment, since I know how to get R-squared and such. The confidence intervals come out really similar to the Bayesian regression.

# Generate dataframes for model fitting
dat_all_gg <- tmgg %>% filter(trait == 'logweight') %>%
  dplyr::select(trait, stat, t, siteID, plotID, bio1, bio4, cv_bio1, bio12, bio15, cv_bio12, chao1, mpd.obs.z2015, mntd.obs.z2015, ruggedness) %>%
  filter(complete.cases(.)) %>%
  rename(mpd_z = mpd.obs.z2015, mntd_z = mntd.obs.z2015) %>%
  mutate_each(funs((. - mean(.))/sd(.)), bio1:ruggedness) 

dat_all_gh <- tmgh %>% filter(trait == 'logweight') %>%
  dplyr::select(trait, stat, t, siteID, plotID, bio1, bio4, cv_bio1, bio12, bio15, cv_bio12, chao1, mpd.obs.z2015, mntd.obs.z2015, ruggedness) %>%
  filter(complete.cases(.)) %>%
  rename(mpd_z = mpd.obs.z2015, mntd_z = mntd.obs.z2015) %>%
  mutate_each(funs((. - mean(.))/sd(.)), bio1:ruggedness) 

dat_all_ggh <- tmggh %>% filter(trait == 'logweight') %>%
  dplyr::select(trait, stat, t, siteID, plotID, bio1, bio4, cv_bio1, bio12, bio15, cv_bio12, chao1, mpd.obs.z2015, mntd.obs.z2015, ruggedness) %>%
  filter(complete.cases(.)) %>%
  rename(mpd_z = mpd.obs.z2015, mntd_z = mntd.obs.z2015) %>%
  mutate_each(funs((. - mean(.))/sd(.)), bio1:ruggedness) 

dat_all_gen <- tmgen %>% filter(trait == 'logweight') %>%
  dplyr::select(trait, stat, t, siteID, plotID, bio1, bio4, cv_bio1, bio12, bio15, cv_bio12, chao1, mpd.obs.z2015, mntd.obs.z2015, ruggedness) %>%
  filter(complete.cases(.)) %>%
  rename(mpd_z = mpd.obs.z2015, mntd_z = mntd.obs.z2015) %>%
  mutate_each(funs((. - mean(.))/sd(.)), bio1:ruggedness) 

dat_all_gran <- tmgran %>% filter(trait == 'logweight') %>%
  dplyr::select(trait, stat, t, siteID, plotID, bio1, bio4, cv_bio1, bio12, bio15, cv_bio12, chao1, mpd.obs.z2015, mntd.obs.z2015, ruggedness) %>%
  filter(complete.cases(.)) %>%
  rename(mpd_z = mpd.obs.z2015, mntd_z = mntd.obs.z2015) %>%
  mutate_each(funs((. - mean(.))/sd(.)), bio1:ruggedness) 

dat_all_herb <- tmherb %>% filter(trait == 'logweight') %>%
  dplyr::select(trait, stat, t, siteID, plotID, bio1, bio4, cv_bio1, bio12, bio15, cv_bio12, chao1, mpd.obs.z2015, mntd.obs.z2015, ruggedness) %>%
  filter(complete.cases(.)) %>%
  rename(mpd_z = mpd.obs.z2015, mntd_z = mntd.obs.z2015) %>%
  mutate_each(funs((. - mean(.))/sd(.)), bio1:ruggedness) 

# Run regressions (freakquentist)
library(lme4)
library(MuMIn)

data_list <- list(dat_all_gg, dat_all_gh, dat_all_ggh, dat_all_gen, dat_all_gran, dat_all_herb)
stat_list <- list('T_IP.IC', 'T_IC.IR', 'T_PC.PR')
fit_list <- list()
confint_list <- list()
coeff_list <- list()
rsq_list <- list()

# This takes a couple of minutes to run since there are 18 models to fit and each one needs to run 999 bootstrap iterations to get CIs.
for (dat_i in data_list) {
  for (stat_j in stat_list) {
    fit_list[[length(fit_list) + 1]] <- lmer(t ~ bio1 + chao1 + mpd_z + ruggedness + (1|siteID), 
                                             data = dat_i, subset = stat == stat_j, na.action = 'na.pass', REML=FALSE)
    confint_list[[length(confint_list) + 1]] <- confint(fit_list[[length(fit_list)]], method = 'boot', nsim = 999)
    coeff_list[[length(coeff_list) + 1]] <- summary(fit_list[[length(fit_list)]])$coef
    rsq_list[[length(rsq_list) + 1]] <- r.squaredGLMM(fit_list[[length(fit_list)]])
  }
}

cimins <- sapply(confint_list, '[', 4:7, 1)
cimaxes <- sapply(confint_list, '[', 4:7, 2)
coeffs <- sapply(coeff_list, '[', 2:5, 1)

r2data <- do.call('rbind', rsq_list)
r2data <- data.frame(r2data, community = rep(c('Gen + Gran', 'Gen + Herb', 'Gen + Gran + Herb', 'Gen only', 'Gran only', 'Herb only'), each = 3), stat = c('T_IP.IC','T_IC.IR','T_PC.PR'))

library(reshape2)
cimins <- melt(cimins, varnames = c('predictor','model'), value.name = 'ci_min')
cimaxes <- melt(cimaxes, varnames = c('predictor','model'), value.name = 'ci_max')
coeffs <- melt(coeffs, varnames = c('predictor','model'), value.name = 'coefficient')

fits_dat <- coeffs %>% left_join(cimins) %>% left_join(cimaxes) %>%
  mutate(community = rep(c('Gen + Gran', 'Gen + Herb', 'Gen + Gran + Herb', 'Gen only', 'Gran only', 'Herb only'), each = 12),
         stat = rep(rep(c('T_IP.IC','T_IC.IR','T_PC.PR'), each = 4), times=6))

# Plots of regression coefficients ----------------------------------------

# Change level order for plotting
fits_dat <- left_join(fits_dat, r2data) %>% 
  mutate(predictor = factor(predictor, levels = c('bio1','ruggedness','mpd_z','chao1'), labels = c('Mean annual\ntemperature','Topographic\nruggedness','Phylogenetic\ndiversity','Species\nrichness')),
         stat = factor(stat, levels = c('T_IP.IC','T_IC.IR','T_PC.PR')))

#labels = c('T_IP/IC\n(local)', 'T_IC/IR\n(regional)', 'T_PC/PR\n(regional ignoring ITV)' )

ggplot(filter(fits_dat, !community %in% c('Gran only', 'Herb only')), aes(x = predictor, y = coefficient, ymin = ci_min, ymax = ci_max)) +
  facet_grid(community ~ stat, labeller = labeller(stat = c('T_IP.IC' = 'T(IP/IC)\nlocal', 'T_IC.IR'='T(IC/IR)\nregional', 'T_PC.PR' = 'T(PC/PR)\nregional without itv'))) + 
  geom_hline(yintercept = 0, color = 'blue', linetype = 'dotted') +
  geom_errorbar(width = 0.1) +
  geom_point() +
  geom_text(aes(x = 4.2, y = -1, label = round(R2m,2))) +
  coord_flip() +
  theme_bw() + theme(panel.grid = element_blank())

ggsave('figs/png/tstat_guild_coefficients_all6.png', height=15, width=9, dpi=400)


# Just show generalist+granivore combo 
coefplotgg <- ggplot(filter(fits_dat, community == 'Gen + Gran'), aes(x = predictor, y = coefficient, ymin = ci_min, ymax = ci_max)) +
  facet_grid(. ~ stat, labeller = labeller(stat = c('T_IP.IC' = 'T(IP/IC)\nlocal', 'T_IC.IR'='T(IC/IR)\nregional', 'T_PC.PR' = 'T(PC/PR)\nregional without itv'))) + 
  geom_hline(yintercept = 0, color = 'blue', linetype = 'dotted') +
  geom_errorbar(width = 0.1) +
  geom_point() +
  geom_text(aes(x = 4.2, y = -1, label = round(R2m,2))) +
  coord_flip() +
  theme_bw() + theme(panel.grid = element_blank())

ggsave('figs/png/tstat_genandgrancoefficients.png', coefplotgg, height = 4, width = 8, dpi = 400)



# Regressions with raw values (beta) --------------------------------------

# Run regressions (freakquentist)
library(glmmADMB)
library(betareg)

data_list <- list(dat_all_gg, dat_all_gh, dat_all_ggh, dat_all_gen, dat_all_gran, dat_all_herb)
stat_list <- list('T_IP.IC', 'T_IC.IR', 'T_PC.PR')
fit_list <- list()
confint_list <- list()
coeff_list <- list()
#rsq_list <- list()

for (dat_i in data_list) {
  for (stat_j in stat_list) {
    fit_list[[length(fit_list) + 1]] <- glmmadmb(tstat ~ bio1 + chao1 + mpd_z + ruggedness + (1|siteID), family = 'beta',
                                             data = subset(dat_i, stat == stat_j) %>% transform(siteID=factor(siteID), tstat=t))
    confint_list[[length(confint_list) + 1]] <- confint(fit_list[[length(fit_list)]])
    coeff_list[[length(coeff_list) + 1]] <- summary(fit_list[[length(fit_list)]])$coef
    #rsq_list[[length(rsq_list) + 1]] <- r.squaredGLMM(fit_list[[length(fit_list)]])
  }
}

cimins <- sapply(confint_list, '[', 4:7, 1)
cimaxes <- sapply(confint_list, '[', 4:7, 2)
coeffs <- sapply(coeff_list, '[', 2:5, 1)

r2data <- do.call('rbind', rsq_list)
r2data <- data.frame(r2data, community = rep(c('Gen + Gran', 'Gen + Herb', 'Gen + Gran + Herb', 'Gen only', 'Gran only', 'Herb only'), each = 3), stat = c('T_IP.IC','T_IC.IR','T_PC.PR'))

library(reshape2)
cimins <- melt(cimins, varnames = c('predictor','model'), value.name = 'ci_min')
cimaxes <- melt(cimaxes, varnames = c('predictor','model'), value.name = 'ci_max')
coeffs <- melt(coeffs, varnames = c('predictor','model'), value.name = 'coefficient')

fits_dat <- coeffs %>% left_join(cimins) %>% left_join(cimaxes) %>%
  mutate(community = rep(c('Gen + Gran', 'Gen + Herb', 'Gen + Gran + Herb', 'Gen only', 'Gran only', 'Herb only'), each = 12),
         stat = rep(rep(c('T_IP.IC','T_IC.IR','T_PC.PR'), each = 4), times=6))

# Plots of regression coefficients ----------------------------------------

# Change level order for plotting
fits_dat <- left_join(fits_dat, r2data) %>% 
  mutate(predictor = factor(predictor, levels = c('bio1','ruggedness','mpd_z','chao1'), labels = c('Mean annual\ntemperature','Topographic\nruggedness','Phylogenetic\ndiversity','Species\nrichness')),
         stat = factor(stat, levels = c('T_IP.IC','T_IC.IR','T_PC.PR')))

#labels = c('T_IP/IC\n(local)', 'T_IC/IR\n(regional)', 'T_PC/PR\n(regional ignoring ITV)' )

ggplot(filter(fits_dat, !community %in% c('Gran only', 'Herb only')), aes(x = predictor, y = coefficient, ymin = ci_min, ymax = ci_max)) +
  facet_grid(community ~ stat, labeller = labeller(stat = c('T_IP.IC' = 'T(IP/IC)\nlocal', 'T_IC.IR'='T(IC/IR)\nregional', 'T_PC.PR' = 'T(PC/PR)\nregional without itv'))) + 
  geom_hline(yintercept = 0, color = 'blue', linetype = 'dotted') +
  geom_errorbar(width = 0.1) +
  geom_point() +
  geom_text(aes(x = 4.2, y = -1, label = round(R2m,2))) +
  coord_flip() +
  theme_bw() + theme(panel.grid = element_blank())

ggsave('figs/png/tstat_guild_coefficients_all6.png', height=15, width=9, dpi=400)


# Just show generalist+granivore combo 
coefplotgg <- ggplot(filter(fits_dat, community == 'Gen + Gran'), aes(x = predictor, y = coefficient, ymin = ci_min, ymax = ci_max)) +
  facet_grid(. ~ stat, labeller = labeller(stat = c('T_IP.IC' = 'T(IP/IC)\nlocal', 'T_IC.IR'='T(IC/IR)\nregional', 'T_PC.PR' = 'T(PC/PR)\nregional without itv'))) + 
  geom_hline(yintercept = 0, color = 'blue', linetype = 'dotted') +
  geom_errorbar(width = 0.1) +
  geom_point() +
  geom_text(aes(x = 4.2, y = -1, label = round(R2m,2))) +
  coord_flip() +
  theme_bw() + theme(panel.grid = element_blank())

ggsave('figs/png/tstat_genandgrancoefficients.png', coefplotgg, height = 4, width = 8, dpi = 400)
