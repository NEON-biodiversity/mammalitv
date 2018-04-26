# Run O-stats on beetle length data if possible
# Author: QDR
# Project: NEON ITV
# Created: 05 Oct 2016
# Last modified: 17 Nov 2016

# Modified 20 Apr 2018: Specify paths

# Path to code is GitHub directory
# Path to data is NEON directory on server

data_path <- '/mnt/research/neon'

# Modified 17 Nov: add new length data and update with most recent O-stat calculation.

beetleLength <- read.csv('exported_data/allbeetlelengths.csv', stringsAsFactors = FALSE)
load('allorganismal_latest.r')
load('allsiteplot_latest.r')

library(ggplot2)
library(dplyr)

source('code/analysis/Ostats.r')
source('code/analysis/densityoverlap.r')
source('code/analysis/ostat2longform.r')
source('code/bioclimnames.r')

# Get rid of unidentified individuals
beetleLength <- filter(beetleLength, species_name != ' ') # Gets rid of 300 individuals.

# Define all the beetles in the US as the regional species pool and run O-stats.
beetleSites <- factor(beetleLength$siteID)
beetlelengthmat <- matrix(beetleLength$length_final, ncol=1)
dimnames(beetlelengthmat) <- list(NULL, 'length')

beetle_reg_lengths <- replicate(nlevels(beetleSites), beetlelengthmat, simplify = FALSE)
beetle_reg_spp <- replicate(nlevels(beetleSites), matrix(factor(beetleLength$species_name), ncol = 1), simplify = FALSE)


beetle_Ostats <- with(beetleLength, Ostats(traits = beetlelengthmat, plots = factor(siteID), sp = factor(species_name), reg_pool_traits = beetle_reg_lengths, reg_pool_sp = beetle_reg_spp, nperm = 99, stat = 'weighted.median'))


# merge with the environmental data
beetleo <- ostat2longform(beetle_Ostats)

spatialstat <- read.csv(file.path(data_path, 'external_data/final_external_data/NEON_spatial_stats_highres.csv'))
spstatmeans <- neonplotdata %>% cbind(spatialstat) %>% group_by(siteID) %>% summarize(ruggedness=mean(tri, na.rm=TRUE))

heterodf <- read.csv(file.path(data_path, 'MS1_RodentOverlap/R_data/heterogeneity.csv'), stringsAsFactors = FALSE)


chao <- function(x) {
  xcomm <- table(x$species_name)
  S_obs <- length(xcomm)
  f1 <- sum(xcomm == 1)
  f2 <- sum(xcomm == 2)
  return(data.frame(chao1 = S_obs + (f1 * (f1 - 1)) / (2 * (f2 + 1))))
}

beetlechao1site <- beetleLength %>% group_by(siteID) %>% do(chao(.))

beetleo <- beetleo %>% rename(siteID=site) %>% full_join(neonsitedata) %>% full_join(spstatmeans) %>% full_join(beetlechao1site) %>% full_join(heterodf) %>% filter(chao1 > 1)


# Plot Beetle O-stats -----------------------------------------------------

posestemplocal <- ggplot(beetleo, aes(x=bio1)) + 
  geom_segment(aes(y = ostat_norm_localnull_seslower, yend = ostat_norm_localnull_sesupper, xend=bio1), alpha = 0.5, size = 1.5, color = 'skyblue') +
  geom_point(aes(y = ostat_norm_localnull_ses), size = 1.5) +
  labs(y = 'Overlap Statistic Local SES', x = parse(text=bioclimnames[1])) +
  theme_minimal() + theme(panel.border = element_rect(fill='transparent')) +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'midnightblue', size=0.9) +
  qSubtitle('Overlap statistics for ground beetles versus MAT', 'Standardized Effect Sizes (local)')

posestempregional <- ggplot(beetleo, aes(x=bio1)) +
  geom_segment(aes(y = ostat_norm_regnull_seslower, yend = ostat_norm_regnull_sesupper, xend=bio1), alpha = 0.5, size = 1.5, color = 'goldenrod') +
  geom_point(aes(y = ostat_norm_regnull_ses), size = 1.5) +
  labs(y = 'Overlap Statistic Regional SES', x = parse(text=bioclimnames[1])) +
  theme_minimal() + theme(panel.border = element_rect(fill='transparent')) +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'midnightblue', size=0.9) +
  qSubtitle('Overlap statistics for ground beetles versus MAT', 'Standardized Effect Sizes (regional)')

porawtemp <- ggplot(beetleo, aes(x=bio1)) + 
  geom_segment(aes(y = ostat_norm_localnull_lower, yend = ostat_norm_localnull_upper, xend=bio1), alpha = 0.5, size = 1.5, color = 'skyblue') +
  geom_segment(aes(y = ostat_norm_regnull_lower, yend = ostat_norm_regnull_upper, xend=bio1), alpha = 0.5, size = 1.5, color = 'goldenrod') +
  geom_point(aes(y = ostat_norm), size = 1.5) +
  labs(y = 'Overlap Statistic', x = parse(text=bioclimnames[1])) +
  theme_minimal() + theme(panel.border = element_rect(fill='transparent')) +
  #geom_hline(yintercept = 0, linetype = 'dotted', color = 'midnightblue', size=0.9) +
  qSubtitle('Overlap statistics for ground beetles versus MAT', 'Raw values, local and regional nulls')

porawtempnonull <- ggplot(beetleo, aes(x=bio1)) + 
 # geom_segment(aes(y = ostat_norm_localnull_lower, yend = ostat_norm_localnull_upper, xend=bio1), alpha = 0.5, size = 1.5, color = 'skyblue') +
 # geom_segment(aes(y = ostat_norm_regnull_lower, yend = ostat_norm_regnull_upper, xend=bio1), alpha = 0.5, size = 1.5, color = 'goldenrod') +
  geom_point(aes(y = ostat_norm), size = 1.5) +
  labs(y = 'Overlap Statistic', x = parse(text=bioclimnames[1])) +
  scale_y_continuous(limits = c(0,0.25)) +
  theme_minimal() + theme(panel.border = element_rect(fill='transparent')) +
  #geom_hline(yintercept = 0, linetype = 'dotted', color = 'midnightblue', size=0.9) +
  ggtitle('Overlap statistics for ground beetles versus MAT')

# Visualize beetle dist with density plot ---------------------------------

#good_sites <- unique(o2015$siteID)[-(4:5)]
sites_temporder <- neonsitedata %>% arrange(bio1) %>% dplyr::select(siteID, bio1)
theme_john <- theme_bw() + theme(panel.grid = element_blank(), 
                                 axis.text = element_text(size = 12), 
                                 axis.title = element_text(size = 18),
                                 text = element_text(family = 'Helvetica'))


beetleplotdat <- beetleLength %>% filter(siteID %in% neonsitedata$siteID, siteID %in% beetleo$siteID, !siteID %in% c('WOOD','BART','KONZ')) %>% mutate(siteID = factor(siteID, levels=sites_temporder$siteID))

pdensbeetle <- ggplot(beetleplotdat) +
  stat_density(adjust = 2, size = 1, aes(x = length_final, group = species_name), fill = 'black', alpha = 0.25, geom='polygon', position = 'identity') + facet_wrap(~ siteID) +
  scale_x_continuous(name = 'body length (mm)', limits = c(-10, 50)) +
  scale_y_continuous(limits = c(0, 1.2), expand = c(0,0)) +
  geom_text(aes(label = round(ostat_norm,3), x = 25, y = 1), color = 'black', data = beetleo %>% filter(siteID %in% neonsitedata$siteID, !siteID %in% c('WOOD','BART','KONZ')) %>% mutate(siteID = factor(siteID, levels=sites_temporder$siteID)) %>% select(siteID, ostat_norm), family = 'Helvetica') +
  geom_text(aes(label = paste0(round(bio1, 2), '°C'), x = 5, y = 1), color = 'black', data = neonsitedata %>% filter(siteID %in% beetleo$siteID, !siteID %in% c('WOOD','BART','KONZ')) %>% mutate(siteID = factor(siteID, levels=sites_temporder$siteID)) %>% select(siteID, bio1), family = 'Helvetica') +
  theme_john + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())


# Beetle Jitter Plot ------------------------------------------------------

beetleo <- beetleo %>%
  mutate(local_significant = ostat_norm_localnull_ses < ostat_norm_localnull_seslower | ostat_norm_localnull_ses > ostat_norm_localnull_sesupper,
         reg_significant = ostat_norm_regnull_ses < ostat_norm_regnull_seslower | ostat_norm_regnull_ses > ostat_norm_regnull_sesupper)

jitterplotdat <- beetleo %>% 
  select(siteID, ostat_norm_localnull_ses, ostat_norm_regnull_ses, local_significant, reg_significant)

jitterplotdat <- with(jitterplotdat, data.frame(siteID=siteID, 
                                                ses = c(ostat_norm_localnull_ses, ostat_norm_regnull_ses),
                                                significant = c(local_significant, reg_significant),
                                                nullmodel = rep(c('Local','Regional'), each=nrow(jitterplotdat))))


jitterplottext <- data.frame(lab = c('More partitioning\nthan expected', 'Neutral', 'More overlap\nthan expected'),
                             x = c(1.5, 1.5, 1.5),
                             y = c(-10, 1, 10))

pj <- ggplot(jitterplotdat, aes(x=nullmodel,y=ses)) +
  geom_hline(yintercept=0, linetype='dotted', color = 'blue', size=1) +
  geom_jitter(aes(color=significant), height=0, width=0.25) +
  geom_text(aes(x,y,label=lab), data=jitterplottext, family = 'Helvetica') +
  scale_x_discrete(name =  'Null model', labels = c('Local','Regional')) +
  scale_y_continuous(name = expression(paste('SES of NO'[local]))) +
  scale_color_manual(values = c('gray75', 'black')) +
  theme_john + theme(legend.position = c(0.88,0.1))


# Beta regression for beetle o-stats --------------------------------------

library(betareg)

beetleo_regdat <- beetleo %>% select(siteID, ostat_norm, bio1, chao1, pc1_productivityheterogeneity, pc2_topographyheterogeneity) %>%
  filter(!is.na(bio1), !is.na(ostat_norm))

reglocal <- betareg(ostat_norm ~ bio1 + chao1 + pc1_productivityheterogeneity + pc2_topographyheterogeneity, link = 'logit', data = beetleo_regdat, na.action = 'na.pass') # Too many predictors

reglocalbio <- betareg(ostat_norm ~ bio1, link = 'logit', data = beetleo_regdat, na.action = 'na.pass')
reglocalchao <- betareg(ostat_norm ~ chao1, link = 'logit', data = beetleo_regdat, na.action = 'na.pass')
reglocalpc1 <- betareg(ostat_norm ~ pc1_productivityheterogeneity, link = 'logit', data = beetleo_regdat, na.action = 'na.pass')
reglocalpc2 <- betareg(ostat_norm ~ pc2_topographyheterogeneity, link = 'logit', data = beetleo_regdat, na.action = 'na.pass')
