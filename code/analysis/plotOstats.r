# Plot O-stats
# Author: QDR
# Project: NEON ITV
# Created: 3 Oct 2016
# Last modified: 31 Oct 2016

# Modified 20 Apr 2018: Specify paths to data.

# Modified 31 Oct: add predictive interval to plot
# Modified 28 Oct: add analysis of residuals for no-itv regression.
# Modified 27 Oct: Add comparison with no-itv case
# Modified 26 Oct: Add the new file names and the updated species richness.
# Modified 19 Oct: Precipitation and seasonality added as predictors in beta reg.
# Modified 18 Oct: Scatterplots with site labels
# Modified 17 Oct: Load O-stats calculated as weighted medians.
# Modified 12 Oct: Replaced the O-stats objects calculated with all individuals, with those calculated for just generalists.
# Modified 10 Oct: Added the regional stats calculated with an entire-continent regional species pool for each site.
# Modified 6 Oct: Added the longer runs of the null models (999 iterations instead of 99) for more accurate confidence intervals. Also the regional o-stats.

data_path <- '/mnt/research/neon'

load(file.path(data_path, 'MS1_RodentOverlap/R_data/overlapstat/Ostats_bysite2014.r'))
load(file.path(data_path, 'MS1_RodentOverlap/R_data/overlapstat/Ostats_bysite2015.r'))
load(file.path(data_path, 'MS1_RodentOverlap/R_data/overlapstat/Ostats_bysiteallyears.r'))

load(file.path(data_path, 'MS1_RodentOverlap/R_data/overlapstat/Ostats_allregional.r'))

# Generalists: load files
files2load <- dir(file.path(data_path, 'MS1_RodentOverlap/R_data/overlapstat'), pattern='gen*')
for (i in file.path(data_path, 'MS1_RodentOverlap/R_data/overlapstat', files2load)) load(i)

# Weighted medians: load files
files2load <- dir(file.path(data_path, 'MS1_RodentOverlap/R_data/overlapstat'), pattern='*wmed.r')
for (i in file.path(data_path, 'MS1_RodentOverlap/R_data/overlapstat', files2load)) load(i)


source('code/analysis/ostat2longform.r')

o2014 <- ostat2longform(Ostats_bysite2014)
o2015 <- ostat2longform(Ostats_bysite2015)
oallyears <- ostat2longform(Ostats_bysiteallyears)

or2014 <- regionalostat2longform(Ostats_regional2014)
or2015 <- regionalostat2longform(Ostats_regional2015)
orallyears <- regionalostat2longform(Ostats_regionalallyears)

orc2015 <- regionalostat2longform(Ostats_regcontinent2015)

# for generalists
o2015 <- ostat2longform(genOstats_bysite2015)
oallyears <- ostat2longform(genOstats_bysiteallyears)
or2015 <- regionalostat2longform(genOstats_regional2015)
orallyears <- regionalostat2longform(genOstats_regionalallyears)
orc2015 <- regionalostat2longform(genOstats_regcontinent2015)
orcallyears <- regionalostat2longform(genOstats_regcontinentallyears)

library(dplyr)
library(ggplot2)
library(lubridate)
library(reshape2)

source('~/GitHub/NEON/code/bioclimnames.r')

load(file.path(data_path, 'MS1_RodentOverlap/R_data/mammalPDbyplotobject.r'))

source('~/GitHub/NEON/code/data_extraction/loadmammallocal.r')
source('~/GitHub/NEON/code/data_extraction/createregpoollists.r')

spatialstat <- read.csv(file.path(data_path, 'external_data/final_external_data/NEON_spatial_stats_highres.csv'))
spstatmeans <- neonplotdata %>% cbind(spatialstat) %>% group_by(siteID) %>% summarize(ruggedness=mean(tri, na.rm=TRUE))

richnessdf <- read.csv('richness.csv', stringsAsFactors = FALSE)

# chao <- function(x) {
#   xcomm <- table(x$taxonID)
#   S_obs <- length(xcomm)
#   f1 <- sum(xcomm == 1)
#   f2 <- sum(xcomm == 2)
#   return(data.frame(chao1 = S_obs + (f1 * (f1 - 1)) / (2 * (f2 + 1))))
# }
# 
# chao1site <- mam_capture %>% group_by(siteID) %>% do(chao(.))
# chao15site <- mam_capture %>% filter(year==2015) %>% group_by(siteID) %>% do(chao(.))
mammalPDsite <- mammalPD %>% mutate(siteID = substr(plotID,1,4)) %>% group_by(siteID) %>% summarize(mpd_z = mean(mpd.obs.z2015, na.rm=T), mntd_z = mean(mntd.obs.z2015, na.rm=T))

o2014 <- o2014 %>% rename(siteID=site) %>% left_join(neonsitedata) %>% left_join(spstatmeans) %>% left_join(mammalPDsite) %>% left_join(richnessdf) %>% filter(chao1 > 1)
o2015 <- o2015 %>% rename(siteID=site) %>% left_join(neonsitedata) %>% left_join(spstatmeans) %>% left_join(mammalPDsite) %>% left_join(richnessdf) %>% filter(chao1 > 1)
oallyears <- oallyears %>% rename(siteID=site) %>% left_join(neonsitedata) %>% left_join(spstatmeans) %>% left_join(mammalPDsite) %>% left_join(richnessdf) %>% filter(chao1 > 1)

or2014 <- or2014 %>% rename(siteID=site) %>% left_join(neonsitedata) %>% left_join(spstatmeans) %>% left_join(mammalPDsite) %>% left_join(richnessdf) %>% filter(chao1 > 1)
or2015 <- or2015 %>% rename(siteID=site) %>% left_join(neonsitedata) %>% left_join(spstatmeans) %>% left_join(mammalPDsite) %>% left_join(richnessdf) %>% filter(chao1 > 1)
orallyears <- orallyears %>% rename(siteID=site) %>% left_join(neonsitedata) %>% left_join(spstatmeans) %>% left_join(mammalPDsite) %>% left_join(richnessdf) %>% filter(chao1 > 1)

orc2015 <- orc2015 %>% rename(siteID=site) %>% left_join(neonsitedata) %>% left_join(spstatmeans) %>% left_join(mammalPDsite) %>% left_join(richnessdf) %>% filter(chao1 > 1)
orcallyears <- orcallyears %>% rename(siteID=site) %>% left_join(neonsitedata) %>% left_join(spstatmeans) %>% left_join(mammalPDsite) %>% left_join(richnessdf) %>% filter(chao1 > 1)

# Scatterplots with null distributions (local and reg) --------------------

o2015goodsites <- filter(o2015, !siteID %in% c('DELA','DSNY','HEAL'))

posestemplocal <- ggplot(o2015goodsites %>% filter(trait=='logweight'), aes(x=bio1)) + 
  geom_segment(aes(y = ostat_norm_localnull_seslower, yend = ostat_norm_localnull_sesupper, xend=bio1), alpha = 0.5, size = 1.5, color = 'skyblue') +
  geom_point(aes(y = ostat_norm_localnull_ses), size = 1.5) +
  labs(y = 'Overlap Statistic Local SES', x = parse(text=bioclimnames[1])) +
  theme_minimal() + theme(panel.border = element_rect(fill='transparent')) +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'midnightblue', size=0.9) +
  qSubtitle('Overlap statistics for 2015 NEON mammals versus MAT', 'Standardized Effect Sizes (local)')

posestempregional <- ggplot(o2015goodsites %>% filter(trait=='logweight'), aes(x=bio1)) +
  geom_segment(aes(y = ostat_norm_regnull_seslower, yend = ostat_norm_regnull_sesupper, xend=bio1), alpha = 0.5, size = 1.5, color = 'goldenrod') +
  geom_point(aes(y = ostat_norm_regnull_ses), size = 1.5) +
  labs(y = 'Overlap Statistic Regional SES', x = parse(text=bioclimnames[1])) +
  theme_minimal() + theme(panel.border = element_rect(fill='transparent')) +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'midnightblue', size=0.9) +
  qSubtitle('Overlap statistics for 2015 NEON mammals versus MAT', 'Standardized Effect Sizes (regional)')

porawtemp <- ggplot(o2015goodsites %>% filter(trait=='logweight'), aes(x=bio1)) + 
  geom_segment(aes(y = ostat_norm_localnull_lower, yend = ostat_norm_localnull_upper, xend=bio1), alpha = 0.5, size = 1.5, color = 'skyblue') +
  geom_segment(aes(y = ostat_norm_regnull_lower, yend = ostat_norm_regnull_upper, xend=bio1), alpha = 0.5, size = 1.5, color = 'goldenrod') +
  geom_point(aes(y = ostat_norm), size = 1.5) +
  labs(y = 'Overlap Statistic', x = parse(text=bioclimnames[1])) +
  theme_minimal() + theme(panel.border = element_rect(fill='transparent')) +
  #geom_hline(yintercept = 0, linetype = 'dotted', color = 'midnightblue', size=0.9) +
  qSubtitle('Overlap statistics for 2015 NEON mammals versus MAT', 'Raw values, local and regional nulls')

poseschaolocal <- ggplot(o2015goodsites %>% filter(trait=='logweight'), aes(x=chao1)) + 
  geom_segment(aes(y = ostat_norm_localnull_seslower, yend = ostat_norm_localnull_sesupper, xend=chao1), alpha = 0.5, size = 1.5, color = 'skyblue') +
  geom_point(aes(y = ostat_norm_localnull_ses), size = 1.5) +
  labs(y = 'Overlap Statistic Local SES', x = 'Species Richness (Chao1)') +
  theme_minimal() + theme(panel.border = element_rect(fill='transparent')) +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'midnightblue', size=0.9) +
  qSubtitle('Overlap statistics for 2015 NEON mammals versus Richness', 'Standardized Effect Sizes (local)')

poseschaoregional <- ggplot(o2015goodsites %>% filter(trait=='logweight'), aes(x=chao1)) +
  geom_segment(aes(y = ostat_norm_regnull_seslower, yend = ostat_norm_regnull_sesupper, xend=chao1), alpha = 0.5, size = 1.5, color = 'goldenrod') +
  geom_point(aes(y = ostat_norm_regnull_ses), size = 1.5) +
  labs(y = 'Overlap Statistic Regional SES', x = 'Species Richness (Chao1)') +
  theme_minimal() + theme(panel.border = element_rect(fill='transparent')) +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'midnightblue', size=0.9) +
  qSubtitle('Overlap statistics for 2015 NEON mammals versus Richness', 'Standardized Effect Sizes (regional)')

porawchao <- ggplot(o2015goodsites %>% filter(trait=='logweight'), aes(x=chao1)) + 
  geom_segment(aes(y = ostat_norm_localnull_lower, yend = ostat_norm_localnull_upper, xend=chao1), alpha = 0.5, size = 1.5, color = 'skyblue') +
  geom_segment(aes(y = ostat_norm_regnull_lower, yend = ostat_norm_regnull_upper, xend=chao1), alpha = 0.5, size = 1.5, color = 'goldenrod') +
  geom_point(aes(y = ostat_norm), size = 1.5) +
  labs(y = 'Overlap Statistic', x = 'Species Richness (Chao1)') +
  theme_minimal() + theme(panel.border = element_rect(fill='transparent')) +
  #geom_hline(yintercept = 0, linetype = 'dotted', color = 'midnightblue', size=0.9) +
  qSubtitle('Overlap statistics for 2015 NEON mammals versus Richness', 'Raw values, local and regional nulls')

posesmpdlocal <- ggplot(o2015goodsites %>% filter(trait=='logweight'), aes(x=mpd_z)) + 
  geom_segment(aes(y = ostat_norm_localnull_seslower, yend = ostat_norm_localnull_sesupper, xend=mpd_z), alpha = 0.5, size = 1.5, color = 'skyblue') +
  geom_point(aes(y = ostat_norm_localnull_ses), size = 1.5) +
  labs(y = 'Overlap Statistic Local SES', x = 'Mean Pairwise Distance SES') +
  theme_minimal() + theme(panel.border = element_rect(fill='transparent')) +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'midnightblue', size=0.9) +
  qSubtitle('Overlap statistics for 2015 NEON mammals versus MPD', 'Standardized Effect Sizes (local)')

posesmpdregional <- ggplot(o2015goodsites %>% filter(trait=='logweight'), aes(x=mpd_z)) +
  geom_segment(aes(y = ostat_norm_regnull_seslower, yend = ostat_norm_regnull_sesupper, xend=mpd_z), alpha = 0.5, size = 1.5, color = 'goldenrod') +
  geom_point(aes(y = ostat_norm_regnull_ses), size = 1.5) +
  labs(y = 'Overlap Statistic Regional SES', x = 'Mean Pairwise Distance SES') +
  theme_minimal() + theme(panel.border = element_rect(fill='transparent')) +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'midnightblue', size=0.9) +
  qSubtitle('Overlap statistics for 2015 NEON mammals versus MPD', 'Standardized Effect Sizes (regional)')

porawmpd <- ggplot(o2015goodsites %>% filter(trait=='logweight'), aes(x=mpd_z)) + 
  geom_segment(aes(y = ostat_norm_localnull_lower, yend = ostat_norm_localnull_upper, xend=mpd_z), alpha = 0.5, size = 1.5, color = 'skyblue') +
  geom_segment(aes(y = ostat_norm_regnull_lower, yend = ostat_norm_regnull_upper, xend=mpd_z), alpha = 0.5, size = 1.5, color = 'goldenrod') +
  geom_point(aes(y = ostat_norm), size = 1.5) +
  labs(y = 'Overlap Statistic', x = 'Mean Pairwise Distance SES') +
  theme_minimal() + theme(panel.border = element_rect(fill='transparent')) +
  #geom_hline(yintercept = 0, linetype = 'dotted', color = 'midnightblue', size=0.9) +
  qSubtitle('Overlap statistics for 2015 NEON mammals versus MPD', 'Raw values, local and regional nulls')


pdf('figs/pdf/ostats_with_nullmodels.pdf', height=6, width=6)
  posestemplocal
  posestempregional
  porawtemp
  poseschaolocal
  poseschaoregional
  porawchao
  posesmpdlocal
  posesmpdregional
  porawmpd
dev.off()

# Scatterplot with site labels
library(directlabels)
sitelab <- geom_dl(aes(label = siteID, y = ostat_norm, size = 5), method = list('top.bumptwice', cex = 0.75, vjust = -0.5))
ggsave('figs/png/ostatscatter/localtemp.png', porawtemp + sitelab, height=6, width=6, dpi=400)
ggsave('figs/png/ostatscatter/localrichness.png', porawchao + sitelab, height=6, width=6, dpi=400)
ggsave('figs/png/ostatscatter/localmpd.png', porawmpd + sitelab, height=6, width=6, dpi=400)

# Determine percent generalist of each community

mam_guildtables <- mam_capture_sitemerge %>% 
  group_by(siteID) %>%
  do(guildtable = table(.$Pineda_Main_food, exclude = 'Carnivore'))

mam_guildtables <- with(mam_guildtables, data.frame(siteID=siteID, do.call('rbind', guildtable))) %>%
  group_by(siteID) %>%
  summarize(pctgen = Generalist / (Generalist+Granivore+Herbivore+Insectivore))

porawgen <- ggplot(o2015goodsites %>% filter(trait=='logweight') %>% left_join(mam_guildtables), aes(x=pctgen)) + 
  geom_segment(aes(y = ostat_norm_localnull_lower, yend = ostat_norm_localnull_upper, xend=pctgen), alpha = 0.5, size = 1.5, color = 'skyblue') +
  geom_segment(aes(y = ostat_norm_regnull_lower, yend = ostat_norm_regnull_upper, xend=pctgen), alpha = 0.5, size = 1.5, color = 'goldenrod') +
  geom_point(aes(y = ostat_norm), size = 1.5) +
  labs(y = 'Overlap Statistic', x = 'Proportion Generalists') +
  theme_minimal() + theme(panel.border = element_rect(fill='transparent')) +
  #geom_hline(yintercept = 0, linetype = 'dotted', color = 'midnightblue', size=0.9) +
  qSubtitle('Overlap statistics versus % generalists', 'Raw values, local and regional nulls')

# Biome of each community
neonplot <- read.csv(file.path(data_path, 'raw_data/spatial_data/spatial_data_csvs/TOSplotSpatialData.csv'))
nlcdtables <- neonplot %>% 
  filter(subtype=='mammalGrid') %>%
  group_by(siteID) %>%
  do(nlcdtable = table(.$nlcdClass))

nlcddf <- data.frame(siteID = nlcdtables$siteID, nlcd = sapply(nlcdtables$nlcdtable, function(x) names(x)[which.max(x)]))

nlcddf$biomegrp <- 'woody'
nlcddf$biomegrp[nlcddf$nlcd %in% c('grasslandHerbaceous', 'shrubScrub','cultivatedCrops')] <- 'herbaceous'

jitterpobiome <- ggplot(o2015goodsites %>% filter(trait=='logweight') %>% left_join(nlcddf), aes(x=nlcd)) + 
  geom_jitter(aes(y = ostat_norm), height=0, width=0.3) +
  labs(y = 'Overlap Statistic', x = 'Biome') +
  theme_minimal() + theme(panel.border = element_rect(fill='transparent')) +
  #geom_hline(yintercept = 0, linetype = 'dotted', color = 'midnightblue', size=0.9) +
  qSubtitle('Overlap statistics by biome', 'Raw values')

jitterpobiomegrp <- ggplot(o2015goodsites %>% filter(trait=='logweight') %>% left_join(nlcddf), aes(x=biomegrp)) + 
  geom_boxplot(aes(y = ostat_norm)) +
  labs(y = 'Overlap Statistic', x = 'Biome Group') +
  theme_minimal() + theme(panel.border = element_rect(fill='transparent')) +
  #geom_hline(yintercept = 0, linetype = 'dotted', color = 'midnightblue', size=0.9) +
  qSubtitle('Overlap statistics by biome', 'Raw values')

# Effect size distributions -----------------------------------------------

library(reshape2)

# Find "significance"

o2015goodsites <- o2015goodsites %>%
  mutate(local_significant = ostat_norm_localnull_ses < ostat_norm_localnull_seslower | ostat_norm_localnull_ses > ostat_norm_localnull_sesupper,
         reg_significant = ostat_norm_regnull_ses < ostat_norm_regnull_seslower | ostat_norm_regnull_ses > ostat_norm_regnull_sesupper)

jitterplotdat <- o2015goodsites %>% 
  filter(trait == 'logweight') %>%
  select(siteID, ostat_norm_localnull_ses, ostat_norm_regnull_ses, local_significant, reg_significant)

jitterplotdat <- with(jitterplotdat, data.frame(siteID=siteID, 
                                                ses = c(ostat_norm_localnull_ses, ostat_norm_regnull_ses),
                                                significant = c(local_significant, reg_significant),
                                                nullmodel = rep(c('Local','Regional'), each=nrow(jitterplotdat))))

jitterplottext <- data.frame(lab = c('More partitioning\nthan expected', 'Neutral', 'More overlap\nthan expected'),
                             x = c(1.5, 1.5, 1.5),
                             y = c(-10, 1, 10))

ggplot(jitterplotdat, aes(x=nullmodel,y=ses)) +
  geom_hline(yintercept=0, linetype='dotted', color = 'blue', size=1) +
  geom_jitter(aes(color=significant), height=0, width=0.25) +
  geom_text(aes(x,y,label=lab), data=jitterplottext) +
  scale_x_discrete(name =  'Null model', labels = c('Local','Regional')) +
  scale_y_continuous(name = 'SES of Overlap Statistic') +
  scale_color_manual(values = c('gray75', 'black')) +
  theme_minimal() + theme(panel.border = element_rect(fill='transparent'), legend.position = c(0.9,0.1))
 
ggsave('figs/pdf/ostat_jitterplot.pdf', height=5, width=5) 


# Plot the shaded density plots (corrected) -------------------------------

good_sites <- unique(o2015$siteID)[-(4:5)]

pdensshade <- ggplot(filter(mam_capture_sitemerge, year == 2015, !siteID %in% c('HEAL','DELA','DSNY')) %>% mutate(siteID = factor(siteID, levels=sites_temporder$siteID))) +
  stat_density(adjust = 2, size = 1, aes(x = log10(weight), group = taxonID), fill = 'black', alpha = 0.25, geom='polygon', position = 'identity') + facet_wrap(~ siteID) +
  scale_x_continuous(name = expression(paste('log'[10],' body mass')), breaks = c(1, 2), labels = c(10, 100)) +
  geom_text(aes(label = round(ostat_norm,3), x = 2.5, y = 15), color = 'black', data = o2015 %>% filter(!siteID %in% c('HEAL','DELA','DSNY'), trait %in% 'logweight') %>% mutate(siteID = factor(siteID, levels=sites_temporder$siteID))) +
  geom_text(aes(label = paste0(round(bio1, 2), '°C'), x = 1, y = 15), color = 'black', data = neonsitedata %>% filter(siteID %in% good_sites) %>% mutate(siteID = factor(siteID, levels=sites_temporder$siteID))) +
    #geom_text(aes(label = paste(round(bio1, 2), 'C'), x = 1, y = 15), color = 'black', data = neonsitedata %>% filter(siteID %in% good_sites) %>% mutate(siteID = factor(siteID, levels=sites_temporder$siteID))) +
  theme_john + ggtitle('Sites ordered by temperature') + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())

source('code/facetadjust.r')

png('figs/png/tstats/densityplotsbytemperature_shaded.png', height = 10, width = 12, res=400, units='in')
facetAdjust(pdensshade)
dev.off()


# Plot the regional species pools -----------------------------------------

regpoolplotdat <- list()

for (i in 1:length(siteregpoolsp_mam15_iucn)) {
  regpoolplotdat[[i]] <- data.frame(siteID = names(siteregpoollist_mam15_iucn)[i],
                                    taxonID = siteregpoolsp_mam15_iucn[[i]],
                                    logweight = log10(siteregpoollist_mam15_iucn[[i]]$weight))
}
regpoolplotdat <- do.call('rbind',regpoolplotdat)


sites_temporder <- neonsitedata %>% arrange(bio1) %>% dplyr::select(siteID, bio1)

theme_john <- theme_bw() + theme(panel.grid = element_blank(), axis.text = element_text(size = 12), axis.title = element_text(size = 18))

# Find significantly positive and negative for regional stats

o2015regstat <- o2015goodsites %>%
  filter(trait == 'logweight') %>%
  arrange(bio1) %>%
  mutate(reg_significant_partition = ostat_norm_regnull_ses < ostat_norm_regnull_seslower,
         reg_significant_overlap =  ostat_norm_regnull_ses > ostat_norm_regnull_sesupper) %>%
  select(siteID, reg_significant_partition, reg_significant_overlap) 

stattext <- rep('neutral', nrow(o2015regstat))
stattext[o2015regstat$reg_significant_partition] <- 'significantly partitioned'
stattext[o2015regstat$reg_significant_overlap] <- 'significantly overlapping'

o2015regstat$stattext <- stattext

# With significance labels
ggplot(filter(mam_capture_sitemerge, year == 2015, !siteID %in% c('HEAL','DELA','DSNY')) %>% mutate(siteID = factor(siteID, levels=sites_temporder$siteID), logweight=log10(weight))) +
  stat_density(adjust = 2, size = 1, aes(x = logweight, group = taxonID), fill = 'black', alpha = 1, geom = 'polygon', position = 'identity', data = regpoolplotdat %>% filter(!siteID %in% c('HEAL','DELA','DSNY'))) +
  stat_density(adjust = 2, size = 1, aes(x = logweight, group = taxonID), fill = 'skyblue', alpha = 0.7, geom='polygon', position = 'identity') + facet_wrap(~ siteID) +
  geom_text(aes(label = stattext, x = 1.5, y = 15), data = o2015regstat %>% mutate(siteID = factor(siteID, levels=sites_temporder$siteID))) +
  scale_x_continuous(breaks = c(1, 2), labels = c(10, 100)) +
  theme_john + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  qSubtitle('Regional species pools (black) and local communities (blue)', 'Result of regional null model shown with text, sites ordered by temp')
  
ggsave('figs/pdf/regionalspeciespoolsdensity.pdf', height=10, width=12)

# With regional overlap proportion listed
# Plot the proportion of the local site that overlaps with the regional species pool

ggplot(filter(mam_capture_sitemerge, year == 2015, !siteID %in% c('HEAL','DELA','DSNY')) %>% mutate(siteID = factor(siteID, levels=sites_temporder$siteID), logweight=log10(weight))) +
  stat_density(adjust = 2, size = 1, aes(x = logweight, group = taxonID), fill = 'black', alpha = 1, geom = 'polygon', position = 'identity', data = regpoolplotdat %>% filter(!siteID %in% c('HEAL','DELA','DSNY'))) +
  stat_density(adjust = 2, size = 1, aes(x = logweight, group = taxonID), fill = 'skyblue', alpha = 1, geom='polygon', position = 'identity') + facet_wrap(~ siteID) +
  geom_text(aes(label = round(regoverlap,3), x = 1.5, y = 15), data = regoverlapdf %>% filter(!siteID %in% c('DELA','HEAL','DSNY')) %>% mutate(siteID = factor(siteID, levels=sites_temporder$siteID))) +
  scale_x_continuous(breaks = c(1, 2), labels = c(10, 100)) +
  theme_john + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  qSubtitle('Regional species pools (black) and local communities (blue)', 'Proportion of traitspace from regional pool shown, sites ordered by temp')
ggsave('figs/pdf/regionalspeciespoolsdensity_withproportions.pdf', height=10, width=12)


pdenslabels <- ggplot(filter(mam_capture_sitemerge, year == 2015, !siteID %in% c('HEAL','DELA','DSNY')) %>% mutate(siteID = factor(siteID, levels=sites_temporder$siteID), logweight=log10(weight))) +
  stat_density(adjust = 1, size = 1, aes(x = logweight), fill = 'black', alpha = 1, geom = 'polygon', position = 'identity', data = regpoolplotdat %>% filter(!siteID %in% c('HEAL','DELA','DSNY'))) +
  stat_density(adjust = 1, size = 1, aes(x = logweight), fill = 'skyblue', alpha = 0.75, geom='polygon', position = 'identity') + facet_wrap(~ siteID) +
  geom_text(aes(label = paste('NM1:',stattextallpool), x = 1.5, y = 5), data = or2015regstat %>% filter(!siteID %in% c('DELA','HEAL','DSNY')) %>% mutate(siteID = factor(siteID, levels=sites_temporder$siteID))) +
  geom_text(aes(label = paste('NM2:',stattextbysp), x = 1.5, y = 4), data = or2015regstat %>% filter(!siteID %in% c('DELA','HEAL','DSNY')) %>% mutate(siteID = factor(siteID, levels=sites_temporder$siteID))) +
  geom_text(aes(label = round(ostat_reg,3), x = 1.5, y = 3), data = or2015regstat %>% filter(!siteID %in% c('DELA','HEAL','DSNY')) %>% mutate(siteID = factor(siteID, levels=sites_temporder$siteID))) +
  scale_x_continuous(breaks = c(1, 2), labels = c(10, 100), name = expression(paste('log'[10], ' body mass'))) +
  theme_john + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  qSubtitle('Regional species pools (black) and local communities (blue)', 'Significance of "regional overlap stat" shown, sites ordered by temp')

source('code/facetadjust.r')

pdf('figs/pdf/regionalspeciespoolsdensity_withlabels.pdf', height=10, width=12)
  facetAdjust(pdenslabels)
dev.off()

# Continental regional species pool
pdensconti <-  ggplot(filter(mam_capture_sitemerge, year == 2015, !siteID %in% c('HEAL','DELA','DSNY')) %>% mutate(siteID = factor(siteID, levels=sites_temporder$siteID), logweight=log10(weight))) +
  stat_density(adjust = 1, size = 1, aes(x = logweight), fill = 'black', alpha = 1, geom = 'polygon', position = 'identity', data = filter(mam_capture_sitemerge,  !siteID %in% c('HEAL','DELA','DSNY')) %>% select(-siteID) %>% mutate(logweight=log10(weight))) +
  stat_density(adjust = 1, size = 1, aes(x = logweight), fill = 'skyblue', alpha = 0.75, geom='polygon', position = 'identity') + facet_wrap(~ siteID) +
  geom_text(aes(label = paste('NM1:',stattextallpool), x = 1.5, y = 5), data = orc2015regstat %>% filter(!siteID %in% c('DELA','HEAL','DSNY')) %>% mutate(siteID = factor(siteID, levels=sites_temporder$siteID))) +
  geom_text(aes(label = paste('NM2:',stattextbysp), x = 1.5, y = 4), data = orc2015regstat %>% filter(!siteID %in% c('DELA','HEAL','DSNY')) %>% mutate(siteID = factor(siteID, levels=sites_temporder$siteID))) +
  geom_text(aes(label = round(ostat_reg,3), x = 1.5, y = 3), data = orc2015regstat %>% filter(!siteID %in% c('DELA','HEAL','DSNY')) %>% mutate(siteID = factor(siteID, levels=sites_temporder$siteID))) +
  scale_x_continuous(breaks = c(1, 2), labels = c(10, 100), name = expression(paste('log'[10], ' body mass'))) +
  theme_john + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  qSubtitle('Continental species pools (black) and local communities (blue)', 'Significance of "regional overlap stat" shown, sites ordered by temp')

pdf('figs/pdf/continentalspeciespoolsdensity_withlabels.pdf', height=10, width=12)
  facetAdjust(pdensconti)
dev.off()

# Plots of regional o-stats -----------------------------------------------

or2015goodsites <- filter(or2015, !siteID %in% c('DELA','DSNY','HEAL'))
or2015goodsites <- or2015goodsites %>%
  mutate(allpool_significant = ostat_reg_allpoolnull_ses < ostat_reg_allpoolnull_seslower | ostat_reg_allpoolnull_ses > ostat_reg_allpoolnull_sesupper,
         bysp_significant = ostat_reg_byspnull_ses < ostat_reg_byspnull_seslower | ostat_reg_byspnull_ses > ostat_reg_byspnull_sesupper)

porrawtemp <- ggplot(or2015goodsites %>% filter(trait=='logweight'), aes(x=bio1)) + 
  geom_segment(aes(y = ostat_reg_allpoolnull_lower, yend = ostat_reg_allpoolnull_upper, xend=bio1), alpha = 0.5, size = 1.5, color = 'seagreen3') +
  geom_segment(aes(y = ostat_reg_byspnull_lower, yend = ostat_reg_byspnull_upper, xend=bio1), alpha = 0.5, size = 1.5, color = 'plum2') +
  geom_point(aes(y = ostat_reg), size = 1.5) +
  labs(y = 'Regional Overlap Statistic', x = parse(text=bioclimnames[1])) +
  theme_minimal() + theme(panel.border = element_rect(fill='transparent')) +
  #geom_hline(yintercept = 0, linetype = 'dotted', color = 'midnightblue', size=0.9) +
  qSubtitle('Regional O-stats for 2015 NEON mammals versus MAT', 'Raw values, all-pool and by-species nulls')

porrawchao <- ggplot(or2015goodsites %>% filter(trait=='logweight'), aes(x=chao1)) + 
  geom_segment(aes(y = ostat_reg_allpoolnull_lower, yend = ostat_reg_allpoolnull_upper, xend=chao1), alpha = 0.5, size = 1.5, color = 'seagreen3') +
  geom_segment(aes(y = ostat_reg_byspnull_lower, yend = ostat_reg_byspnull_upper, xend=chao1), alpha = 0.5, size = 1.5, color = 'plum2') +
  geom_point(aes(y = ostat_reg), size = 1.5) +
  labs(y = 'Regional Overlap Statistic', x = 'Chao1 Species Richness') +
  theme_minimal() + theme(panel.border = element_rect(fill='transparent')) +
  #geom_hline(yintercept = 0, linetype = 'dotted', color = 'midnightblue', size=0.9) +
  qSubtitle('Regional O-stats for 2015 NEON mammals versus richness', 'Raw values, all-pool and by-species nulls')

porrawmpd <- ggplot(or2015goodsites %>% filter(trait=='logweight'), aes(x=mpd_z)) + 
  geom_segment(aes(y = ostat_reg_allpoolnull_lower, yend = ostat_reg_allpoolnull_upper, xend=mpd_z), alpha = 0.5, size = 1.5, color = 'seagreen3') +
  geom_segment(aes(y = ostat_reg_byspnull_lower, yend = ostat_reg_byspnull_upper, xend=mpd_z), alpha = 0.5, size = 1.5, color = 'plum2') +
  geom_point(aes(y = ostat_reg), size = 1.5) +
  labs(y = 'Regional Overlap Statistic', x = 'SES of mean pairwise PD') +
  theme_minimal() + theme(panel.border = element_rect(fill='transparent')) +
  #geom_hline(yintercept = 0, linetype = 'dotted', color = 'midnightblue', size=0.9) +
  qSubtitle('Regional O-stats for 2015 NEON mammals versus mpd', 'Raw values, all-pool and by-species nulls')

pdf('figs/pdf/regionalostats_scatterplots.pdf', height=6, width=6)
  porrawtemp
  porrawchao
  porrawmpd
dev.off()

# Regional with sites labeled
sitelabreg <- geom_dl(aes(label = siteID, y = ostat_reg), method = list('top.bumptwice', cex = 0.75, vjust = -0.5))
ggsave('figs/png/ostatscatter/regionaltemp.png', porrawtemp + sitelabreg, height=6, width=6, dpi=400)
ggsave('figs/png/ostatscatter/regionalrichness.png', porrawchao + sitelabreg, height=6, width=6, dpi=400)
ggsave('figs/png/ostatscatter/regionalmpd.png', porrawmpd + sitelabreg, height=6, width=6, dpi=400)


or2015regstat <- or2015goodsites %>%
  filter(trait == 'logweight') %>%
  arrange(bio1) %>%
  mutate(bysp_significant_filter = ostat_reg_byspnull_ses < ostat_reg_byspnull_seslower,
         bysp_significant_overdisperse =  ostat_reg_byspnull_ses > ostat_reg_byspnull_sesupper,
         allpool_significant_filter = ostat_reg_allpoolnull_ses < ostat_reg_allpoolnull_seslower,
         allpool_significant_overdisperse = ostat_reg_allpoolnull_ses > ostat_reg_allpoolnull_sesupper) %>%
  select(siteID, ostat_reg, bysp_significant_filter, bysp_significant_overdisperse, allpool_significant_filter, allpool_significant_overdisperse) 

stattext <- rep('neutral', nrow(or2015regstat))
stattext[or2015regstat$bysp_significant_filter] <- 'filtered'
stattext[or2015regstat$bysp_significant_overdisperse] <- 'overdispersed'

stattextallpool <- rep('neutral', nrow(or2015regstat))
stattextallpool[or2015regstat$allpool_significant_filter] <- 'filtered'
stattextallpool[or2015regstat$allpool_significant_overdisperse] <- 'overdispersed'

or2015regstat$stattextbysp <- stattext
or2015regstat$stattextallpool <- stattextallpool

# Continental stats
orc2015goodsites <- filter(orc2015, !siteID %in% c('DELA','DSNY','HEAL'))
orc2015goodsites <- orc2015goodsites %>%
  mutate(allpool_significant = ostat_reg_allpoolnull_ses < ostat_reg_allpoolnull_seslower | ostat_reg_allpoolnull_ses > ostat_reg_allpoolnull_sesupper,
         bysp_significant = ostat_reg_byspnull_ses < ostat_reg_byspnull_seslower | ostat_reg_byspnull_ses > ostat_reg_byspnull_sesupper)

orc2015regstat <- orc2015goodsites %>%
  filter(trait == 'logweight') %>%
  arrange(bio1) %>%
  mutate(bysp_significant_filter = ostat_reg_byspnull_ses < ostat_reg_byspnull_seslower,
         bysp_significant_overdisperse =  ostat_reg_byspnull_ses > ostat_reg_byspnull_sesupper,
         allpool_significant_filter = ostat_reg_allpoolnull_ses < ostat_reg_allpoolnull_seslower,
         allpool_significant_overdisperse = ostat_reg_allpoolnull_ses > ostat_reg_allpoolnull_sesupper) %>%
  select(siteID, ostat_reg, bysp_significant_filter, bysp_significant_overdisperse, allpool_significant_filter, allpool_significant_overdisperse) 
stattext <- rep('neutral', nrow(orc2015regstat))
stattext[orc2015regstat$bysp_significant_filter] <- 'filtered'
stattext[orc2015regstat$bysp_significant_overdisperse] <- 'overdispersed'

stattextallpool <- rep('neutral', nrow(orc2015regstat))
stattextallpool[orc2015regstat$allpool_significant_filter] <- 'filtered'
stattextallpool[orc2015regstat$allpool_significant_overdisperse] <- 'overdispersed'

orc2015regstat$stattextbysp <- stattext
orc2015regstat$stattextallpool <- stattextallpool

# Beta Reg for local and regional -----------------------------------------

library(betareg)

# local o-stat
reglocal <- betareg(ostat_norm ~ bio1 + chao1 + mpd_z + ruggedness, link = 'logit', data = o2015 %>% filter(trait == 'logweight', !siteID %in% c('DSNY','DELA','HARV')))

summary(reglocal)
confint(reglocal)

reglocalbio <- betareg(ostat_norm ~ bio1, link = 'logit', data = o2015 %>% filter(trait == 'logweight', !siteID %in% c('DSNY','DELA','HARV')))
reglocalchao <- betareg(ostat_norm ~ chao1, link = 'logit', data = o2015 %>% filter(trait == 'logweight', !siteID %in% c('DSNY','DELA','HARV')))
reglocalmpd <- betareg(ostat_norm ~ mpd_z, link = 'logit', data = o2015 %>% filter(trait == 'logweight', !siteID %in% c('DSNY','DELA','HARV')))

# regional o-stat

regregional <- betareg(ostat_reg ~ bio1 + chao1 + mpd_z + ruggedness, link = 'logit', data = or2015 %>% filter(trait == 'logweight', !siteID %in% c('DSNY','DELA','HARV')))

summary(regregional)
confint(regregional)

regregionalbio <- betareg(ostat_reg ~ bio1, link = 'logit', data = or2015 %>% filter(trait == 'logweight', !siteID %in% c('DSNY','DELA')))
regregionalchao <- betareg(ostat_reg ~ chao1, link = 'logit', data = or2015 %>% filter(trait == 'logweight', !siteID %in% c('DSNY','DELA')))
regregionalmpd <- betareg(ostat_reg ~ mpd_z, link = 'logit', data = or2015 %>% filter(trait == 'logweight', !siteID %in% c('DSNY','DELA')))

# regional o-stat with continental pool
regconti <- betareg(ostat_reg ~ bio1 + chao1 + mpd_z + ruggedness, link = 'logit', data = orc2015 %>% filter(trait == 'logweight', !siteID %in% c('DSNY','DELA')))

summary(regconti)
confint(regconti)

regcontibio <- betareg(ostat_reg ~ bio1, link = 'logit', data = orc2015 %>% filter(trait == 'logweight', !siteID %in% c('DSNY','DELA')))
regcontichao <- betareg(ostat_reg ~ chao1, link = 'logit', data = orc2015 %>% filter(trait == 'logweight', !siteID %in% c('DSNY','DELA')))
regcontimpd <- betareg(ostat_reg ~ mpd_z, link = 'logit', data = orc2015 %>% filter(trait == 'logweight', !siteID %in% c('DSNY','DELA')))

# Beta reg with precipitation, local and regional (added 19 Oct.)
reglocalprec <- betareg(ostat_norm ~ bio1 + bio12 + bio4 + bio15 + cv_bio1 + cv_bio12, link = 'logit', data = o2015 %>% filter(trait == 'logweight', !siteID %in% c('DSNY','DELA')))

summary(reglocalprec)

regregionalprec <- betareg(ostat_reg ~ bio1 + bio12 + bio4 + bio15 + cv_bio1 + cv_bio12, link = 'logit', data = or2015 %>% filter(trait == 'logweight', !siteID %in% c('DSNY','DELA')))

summary(regregionalprec)


# Comparison with noITV
reglocalnoitv <- lm(ostat_norm ~ bio1 + chao1 + mpd_z + ruggedness, data = noitv2015 %>% filter(trait == 'logweight', !siteID %in% c('DSNY','DELA'))) 
summary(reglocalnoitv)
confint(reglocalnoitv)

bidat <- data.frame(localmeandiff=subset(noitv2015, trait=='logweight' & !siteID %in% c('DSNY','DELA'))$ostat_norm, localoverlap=subset(o2015, trait=='logweight' & !siteID %in% c('DSNY','DELA'))$ostat_norm)
ggplot(bidat, aes(x=localmeandiff, y=localoverlap)) + geom_point() + theme_john + labs(y='Local overlap', x='Median pairwise mean difference')
ggsave('figs/msfigs/itv_versus_noitv.png', height=5, width=5, dpi=400)

# 28 Oct: see if variance goes down with mean going up

bidat <- bidat[complete.cases(bidat), ]
bilm <- lm(localoverlap~localmeandiff, data=bidat)
bibeta <- betareg(localoverlap~localmeandiff, data=bidat, link='logit')
plot(bilm)
plot(bidat$localmeandiff, abs(bilm$residuals))
plot(bidat$localmeandiff, abs(bibeta$residuals))
residlm <- lm(abs(bibeta$residuals) ~ bidat$localmeandiff)
ggplot(data.frame(localmeandiff=bidat$localmeandiff, abs_residual=abs(bibeta$residuals)), aes(localmeandiff,abs_residual)) + 
  geom_point() + theme_john + labs(y = 'Absolute value of residual', x = 'Median pairwise mean difference') +
  stat_smooth(method='lm', color='blue',size=1.5, se=F)
ggsave('figs/msfigs/residual_vs_noitv.png', height=5, width=5, dpi=400)

predbeta <- predict(bibeta, newdata=data.frame(localmeandiff=seq(0,1.55,length.out=10001)), type='quantile', at=c(0.025, 0.975))
predbeta <- data.frame(cilow=predbeta[,1], cihigh=predbeta[,2],localmeandiff= seq(0,1.55,length.out=10001))

ggplot(bidat) + geom_point(aes(x=localmeandiff, y=localoverlap), size=3) + 
  stat_function(aes(x=x), data=data.frame(x=c(0,1.55)), geom='line', fun=fx, args=list(b0=.7238, b1=-3.2147), color='blue', size=1.2, n = 10001) +
  geom_line(data=predbeta, aes(y=cilow, x=localmeandiff)) +
  geom_line(data=predbeta, aes(y=cihigh, x=localmeandiff)) +
  theme_john + labs(y=expression(NO[local]), x='Median pairwise mean difference') +
  scale_x_continuous(expand=c(0,0), limits=c(0,1.55), breaks=c(0,0.5,1,1.5)) + scale_y_continuous(limits=c(0,1), expand=c(0.01,0), breaks=c(0, 0.5, 1))
ggsave('figs/msfigs/itv_versus_noitv.png', height=5, width=5, dpi=600)


# 27 Oct: add the pca's to the regression instead of ruggedness
reglocal <- betareg(ostat_norm ~ bio1 + chao1 + mpd_z + pc1_productivityheterogeneity + pc2_topographyheterogeneity, link = 'logit', data = o2015 %>% filter(trait == 'logweight', !siteID %in% c('DSNY','DELA', 'HARV')), na.action = 'na.pass')

library(MuMIn)
dredge(reglocal) # Results in best model including only temperature and richness
reglocal2 <- betareg(ostat_norm ~ bio1 + chao1, link = 'logit', data = o2015 %>% filter(trait == 'logweight', !siteID %in% c('DSNY','DELA', 'HARV')), na.action = 'na.pass')
reglocal2_best4 <- betareg(ostat_norm ~ bio1 + chao1 + mpd_z + pc1_productivityheterogeneity, link = 'logit', data = o2015 %>% filter(trait == 'logweight', !siteID %in% c('DSNY','DELA', 'HARV')), na.action = 'na.pass')


regregional <- betareg(ostat_reg ~ bio1 + chao1 + mpd_z + pc1_productivityheterogeneity + pc2_topographyheterogeneity, link = 'logit', data = or2015 %>% filter(trait == 'logweight', !siteID %in% c('DSNY','DELA')), na.action = 'na.pass')

dredge(regregional) # Results in best model including only richness and phylogenetic distance
regregional2 <- betareg(ostat_reg ~ chao1 + mpd_z, link = 'logit', data = or2015 %>% filter(trait == 'logweight', !siteID %in% c('DSNY','DELA')), na.action = 'na.pass')


# regional o-stat with continental pool
regconti <- betareg(ostat_reg ~ bio1 + chao1 + mpd_z + pc1_productivityheterogeneity + pc2_topographyheterogeneity, link = 'logit', data = orc2015 %>% filter(trait == 'logweight', !siteID %in% c('DSNY','DELA')), na.action = 'na.pass')

summary(regconti)
confint(regconti)

library(MuMIn)
dredge(regconti)

regcontibio <- betareg(ostat_reg ~ bio1, link = 'logit', data = orc2015 %>% filter(trait == 'logweight', !siteID %in% c('DSNY','DELA')))
regcontichao <- betareg(ostat_reg ~ chao1, link = 'logit', data = orc2015 %>% filter(trait == 'logweight', !siteID %in% c('DSNY','DELA')))
regcontimpd <- betareg(ostat_reg ~ mpd_z, link = 'logit', data = orc2015 %>% filter(trait == 'logweight', !siteID %in% c('DSNY','DELA')))

# 06 Dec: Different link functions besides logit:
reglocallogit <- betareg(ostat_norm ~ bio1 + chao1 + mpd_z + pc1_productivityheterogeneity + pc2_topographyheterogeneity, link = 'logit', data = o2015 %>% filter(trait == 'logweight', !siteID %in% c('DSNY','DELA', 'HARV')), na.action = 'na.pass')

reglocalprobit <- betareg(ostat_norm ~ bio1 + chao1 + mpd_z + pc1_productivityheterogeneity + pc2_topographyheterogeneity, link = 'probit', data = o2015 %>% filter(trait == 'logweight', !siteID %in% c('DSNY','DELA', 'HARV')), na.action = 'na.pass')

reglocalclog <- betareg(ostat_norm ~ bio1 + chao1 + mpd_z + pc1_productivityheterogeneity + pc2_topographyheterogeneity, link = 'cloglog', data = o2015 %>% filter(trait == 'logweight', !siteID %in% c('DSNY','DELA', 'HARV')), na.action = 'na.pass')

reglocalcauchit <- betareg(ostat_norm ~ bio1 + chao1 + mpd_z + pc1_productivityheterogeneity + pc2_topographyheterogeneity, link = 'cauchit', data = o2015 %>% filter(trait == 'logweight', !siteID %in% c('DSNY','DELA', 'HARV')), na.action = 'na.pass')

reglocalloglog <- betareg(ostat_norm ~ bio1 + chao1 + mpd_z + pc1_productivityheterogeneity + pc2_topographyheterogeneity, link = 'loglog', data = o2015 %>% filter(trait == 'logweight', !siteID %in% c('DSNY','DELA', 'HARV')), na.action = 'na.pass')

library(MuMIn)
dredge(reglocalloglog)
reglocalloglogbest <- betareg(ostat_norm ~ bio1 + chao1, link = 'loglog', data = o2015 %>% filter(trait == 'logweight', !siteID %in% c('DSNY','DELA', 'HARV')), na.action = 'na.pass')

loglogbio <- betareg(ostat_norm ~ bio1, link = 'loglog', data = o2015 %>% filter(trait == 'logweight', !siteID %in% c('DSNY','DELA', 'HARV')), na.action = 'na.pass')
loglogchao <- betareg(ostat_norm ~ chao1, link = 'loglog', data = o2015 %>% filter(trait == 'logweight', !siteID %in% c('DSNY','DELA', 'HARV')), na.action = 'na.pass')

# diagnostic
hist(reglocalloglogbest$residuals)

reglogittrans <- lm(I(qlogis(ostat_norm)) ~ bio1 + chao1, data = o2015 %>% filter(trait == 'logweight', !siteID %in% c('DSNY','DELA', 'HARV')), na.action = 'na.pass')

# check for multicollinearity

library(usdm)
pred_df <- o2015 %>% filter(trait == 'logweight', !siteID %in% c('DSNY','DELA', 'HARV')) %>% select(bio1, chao1, mpd_z, pc1_productivityheterogeneity, pc2_topographyheterogeneity)
pairs(pred_df)
cor(pred_df)
vif(pred_df)
