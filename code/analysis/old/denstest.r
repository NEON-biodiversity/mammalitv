# density overlap test
# Author: QDR
# Project: NEON ITV
# Created: 27 Sep 2016
# Last modified: 29 Sep 2016

# Modified 20 Apr 2018: Specify path to data

data_path <- '/mnt/research/neon'

source('code/analysis/densityoverlap.r')

# harv <- mam_capture_sitemerge %>% filter(siteID=='HARV', year==2015) %>% select(taxonID, weight) %>% filter(complete.cases(.), taxonID != 'PEME')
# 
# community_overlap(traits=log10(harv$weight), sp=harv$taxonID, norm=T)
# community_overlap(traits=log10(harv$weight), sp=harv$taxonID, norm=F)

# sptables <- mam_capture_sitemerge %>% 
#   filter(year==2015, !is.na(weight)) %>% 
#   mutate(logweight = log10(weight)) %>% 
#   group_by(siteID) %>%
#   do(sptable = table(.$taxonID))

overlapstatsbysite <- mam_capture_sitemerge %>% 
  filter(year==2015, siteID != 'DSNY') %>% # Expunge Disney because of poor sampling there 
  mutate(logweight = log10(weight)) %>% 
  group_by(siteID) %>%
  do(overlap_norm = community_overlap(traits = .$logweight, sp = .$taxonID, norm = TRUE),
            overlap_unnorm = community_overlap(traits = .$logweight, sp = .$taxonID, norm = FALSE)) 

overlapstatsbysite <- with(overlapstatsbysite, data.frame(siteID=siteID, overlap_norm=unlist(overlap_norm), overlap_unnorm=unlist(overlap_unnorm)))

spatialstat <- read.csv(file.path(data_path, 'external_data/final_external_data/NEON_spatial_stats_highres.csv'))
load(file.path(data_path, 'MS1_RodentOverlap/R_data/mammalPDbyplotobject.r'))

spstatmeans <- neonplotdata %>% cbind(spatialstat) %>% group_by(siteID) %>% summarize(ruggedness=mean(tri, na.rm=TRUE))
pdmeans <- mammalPD %>% mutate(siteID = sapply(strsplit(plotID2015,'_'), '[', 1)) %>% group_by(siteID) %>% summarize(mpd_z = mean(mpd.obs.z2015, na.rm=TRUE))

chao <- function(x) {
  xcomm <- table(x$taxonID)
  S_obs <- length(xcomm)
  f1 <- sum(xcomm == 1)
  f2 <- sum(xcomm == 2)
  return(data.frame(chao1 = S_obs + (f1 * (f1 - 1)) / (2 * (f2 + 1))))
}

chao15site <- mam_capture %>% filter(year(date) == 2015) %>%
  group_by(siteID) %>% do(chao(.))

overlapdf <- overlapstatsbysite %>% left_join(neonsitedata) %>% left_join(spstatmeans) %>% left_join(pdmeans) %>% left_join(chao15site)

dat_site_overlap <- overlapdf %>% 
  dplyr::select(overlap_norm, overlap_unnorm, siteID, bio1, bio4, cv_bio1, bio12, bio15, cv_bio12, chao1, mpd_z, ruggedness) %>%
  filter(complete.cases(.)) %>%
  mutate_each(funs((. - mean(.))/sd(.)), bio1:ruggedness) 

plot(overlap_norm ~ bio1, data=overlapdf)
plot(overlap_norm ~ chao1, data=overlapdf)
plot(overlap_norm ~ ruggedness, data=overlapdf)
plot(overlap_norm ~ mpd_z, data=overlapdf)

library(betareg)

reg <- betareg(overlap_norm ~ bio1 + chao1 + mpd_z + ruggedness, link='logit', data=dat_site_overlap)
reg2 <- betareg(overlap_unnorm ~ bio1 + chao1 + mpd_z + ruggedness, link='logit', data=dat_site_overlap)

confint(reg)
confint(reg2)

regbio <- betareg(overlap_norm ~ bio1, link='logit', data=dat_site_overlap)
regchao <- betareg(overlap_norm ~ chao1, link='logit', data=dat_site_overlap)


library(brms)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

bayesianreg <- brm(overlap_norm ~ bio1 + chao1 + mpd_z + ruggedness, family = Beta(link='logit'), data=dat_site_overlap, chains=4, cores=4, iter = 10000, warmup = 5000)
bayesianreg2 <- brm(overlap_unnorm ~ bio1 + chao1 + mpd_z + ruggedness, family = Beta(link='logit'), data=dat_site_overlap, chains=4, cores=4, iter = 10000, warmup = 5000)

summary(bayesianreg)
summary(bayesianreg2)

###############

# Run overlap stats by plot as well.
overlapstatsbyplot <- mam_capture_sitemerge %>% 
  filter(year==2015) %>% 
  mutate(logweight = log10(weight)) %>% 
  group_by(plotID) %>%
  do(overlap_norm = community_overlap(traits = .$logweight, sp = .$taxonID, norm = TRUE),
     overlap_unnorm = community_overlap(traits = .$logweight, sp = .$taxonID, norm = FALSE)) 

overlapstatsbyplot <- with(overlapstatsbyplot, data.frame(plotID=plotID, overlap_norm=unlist(overlap_norm), overlap_unnorm=unlist(overlap_unnorm)))

spatialstat <- read.csv(file.path(data_path, 'external_data/final_external_data/NEON_spatial_stats_highres.csv'))
neonplotdata <- cbind(neonplotdata, spatialstat) %>% rename(ruggedness = tri)
overlapdfbyplot <- left_join(overlapstatsbyplot, filter(neonplotdata, subtype == 'mammalGrid'))

chao <- function(x) {
  xcomm <- table(x$taxonID)
  S_obs <- length(xcomm)
  f1 <- sum(xcomm == 1)
  f2 <- sum(xcomm == 2)
  return(data.frame(chao1 = S_obs + (f1 * (f1 - 1)) / (2 * (f2 + 1))))
}

chao15 <- mam_capture %>% filter(year(date) == 2015) %>%
  group_by(plotID) %>% do(chao(.))

load(file.path(data_path, 'MS1_RodentOverlap/R_data/mammalPDbyplotobject.r'))

overlapdfbyplot <- overlapdfbyplot %>% left_join(chao15) %>% left_join(mammalPD) %>% filter(chao1 > 1)

dat_all_overlap <- overlapdfbyplot %>% 
  dplyr::select(overlap_norm, overlap_unnorm, siteID, plotID, bio1, bio4, cv_bio1, bio12, bio15, cv_bio12, chao1, mpd.obs.z2015, mntd.obs.z2015, ruggedness) %>%
  filter(complete.cases(.)) %>%
  rename(mpd_z = mpd.obs.z2015, mntd_z = mntd.obs.z2015) %>%
  mutate_each(funs((. - mean(.))/sd(.)), bio1:ruggedness) 



rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

byplotreg <- brm(overlap_norm ~ bio1 + mpd_z + chao1 + ruggedness + (1|siteID), family = Beta(link='logit'), data=dat_all_overlap, chains=4, cores=4, iter = 10000, warmup = 5000)
byplotreg2 <- brm(overlap_unnorm ~ bio1 + mpd_z + chao1 + ruggedness + (1|siteID), family = Beta(link='logit'), data=dat_all_overlap, chains=4, cores=4, iter = 10000, warmup = 5000)

summary(byplotreg)
summary(byplotreg2)

library(glmmADMB)

byplotbetareg <- glmmadmb(overlap_norm ~ bio1 + mpd_z + chao1 + ruggedness + (1|siteID), family = 'beta', data=transform(dat_all_overlap, siteID=factor(siteID)))
byplotbetareg2 <- glmmadmb(overlap_unnorm ~ bio1 + mpd_z + chao1 + ruggedness + (1|siteID), family = 'beta', data=transform(dat_all_overlap, siteID=factor(siteID)))

confint(byplotbetareg)
confint(byplotbetareg2)

library(MuMIn)
r.squaredGLMM(byplotbetareg) # Doesn't work.

plot(overlap_norm ~ ruggedness, data=overlapdfbyplot)



# Plots by site -----------------------------------------------------------

library(ggplot2)
library(GGally)

ggpairs(data=overlapdf[,c('overlap_norm','bio1','ruggedness','chao1','mpd_z')])

theme_john <- theme_bw() + theme(panel.grid = element_blank(), axis.text = element_text(size = 12), axis.title = element_text(size = 18))

overlapdf <- filter(overlapdf, !is.na(overlap_norm))

pdf('figs/pdf/overlapstat_by_predictors.pdf', height=5, width=5)

ggplot(overlapdf, aes(bio1, overlap_norm)) + 
  stat_smooth(method='lm',se=F) +
  labs(x = 'Mean Annual Temperature') +
  geom_point() + theme_john
# ggplot(overlapdf, aes(bio12, overlap_norm)) + 
#   stat_smooth(method='lm',se=T) +
#   geom_point() + theme_john
ggplot(overlapdf, aes(ruggedness, overlap_norm)) + 
  labs(x = 'Topographic Ruggedness') +
  geom_point() + theme_john
ggplot(overlapdf, aes(chao1, overlap_norm)) + 
  stat_smooth(method='lm',se=F) +
  labs(x = 'Chao Species Richness') +
  geom_point() + theme_john
ggplot(overlapdf, aes(mpd_z, overlap_norm)) + 
  labs(x = 'Mean Pairwise Phy. Dist. effect size') +
  geom_point() + theme_john

dev.off()



overlapdfbyplot <- filter(overlapdfbyplot, !is.na(overlap_norm))

pdf('figs/pdf/overlapstat_by_predictors_byplot.pdf', height=5, width=5)

ggplot(overlapdfbyplot, aes(bio1, overlap_norm, color=siteID)) + 
  #stat_smooth(method='lm',se=T) +
  labs(x = 'Mean Annual Temperature') +
  geom_point() + theme_john
ggplot(overlapdfbyplot, aes(bio4, overlap_norm)) + 
  stat_smooth(method='lm',se=T) +
  labs(x = 'Temp. Seasonality') +
  geom_point() + theme_john
ggplot(overlapdfbyplot, aes(bio15, overlap_norm)) + 
  stat_smooth(method='lm',se=T) +
  labs(x = 'Precip. Seasonality') +
  geom_point() + theme_john
ggplot(overlapdfbyplot, aes(cv_bio1, overlap_norm)) + 
  stat_smooth(method='lm',se=T) +
  labs(x = 'Temp Interannual Var') +
  geom_point() + theme_john
ggplot(overlapdfbyplot, aes(bio12, overlap_norm)) + 
  stat_smooth(method='lm',se=T) +
  labs(x = 'Yearly Precip') +
  geom_point() + theme_john
ggplot(overlapdfbyplot, aes(ruggedness, overlap_norm)) + 
  stat_smooth(method='lm',se=T) +
  geom_point() + theme_john
ggplot(overlapdfbyplot, aes(chao1, overlap_norm)) + 
  stat_smooth(method='lm',se=T) +
  labs(x = 'Species Richness') +
  geom_point() + theme_john
ggplot(overlapdfbyplot, aes(mpd.obs.z2015, overlap_norm)) + 
  stat_smooth(method='lm',se=T) +
  labs(x = 'Pairwise Phylogenetic Diversity') +
  geom_point() + theme_john

dev.off()

# Overlap distributions by site -------------------------------------------

overlapstatsbysite_dist <- mam_capture_sitemerge %>% 
  filter(year==2015, siteID != 'DSNY') %>% # Expunge Disney because of poor sampling there 
  mutate(logweight = log10(weight)) %>% 
  group_by(siteID) %>%
  do(overlap_norm = community_overlap_distributions(traits = .$logweight, sp = .$taxonID, norm = TRUE),
     overlap_unnorm = community_overlap_distributions(traits = .$logweight, sp = .$taxonID, norm = FALSE)) 

# Plot distributions of pairwise overlaps for each site

# Normalized overlaps
#hist( overlapstatsbysite_dist$overlap_norm[[3]])

# Convert to data frame for plotting
overlapdistdf <- with(overlapstatsbysite_dist, data.frame(siteID = rep(siteID, times = sapply(overlap_norm, length)), 
                                                          overlap_norm = unlist(overlap_norm),
                                                          overlap_unnorm = unlist(overlap_unnorm)))

sites_statorder <- overlapstatsbysite %>% filter(!is.na(overlap_norm)) %>% arrange(overlap_norm)

library(ggplot2)
ggplot(filter(overlapdistdf, !is.na(overlap_norm)) %>% mutate(siteID=factor(siteID, levels=sites_statorder$siteID)), aes(x = overlap_norm)) +
  facet_wrap(~ siteID, scales = 'free_y') +
  geom_histogram(bins=10) + 
  geom_text(aes(label = round(overlap_norm,3)), x = Inf, y = Inf, hjust = 1, vjust = 1, color = 'red', data = sites_statorder) +
  theme_john
ggsave('figs/png/tstats/overlapstatdistributions.png', height=12, width=12, dpi=400)

ggplot(filter(overlapdistdf, !is.na(overlap_norm)) %>% mutate(siteID=factor(siteID, levels=sites_statorder$siteID)), aes(x = overlap_norm)) +
  facet_wrap(~ siteID, scales = 'free_y') +
  geom_density(fill = 'gray50') + 
  geom_text(aes(label = round(overlap_norm,3)), x = Inf, y = Inf, hjust = 1, vjust = 1, color = 'red', data = sites_statorder) +
  theme_john
ggsave('figs/png/tstats/overlapstatdistributions_densityplots.png', height=12, width=12, dpi=400)
