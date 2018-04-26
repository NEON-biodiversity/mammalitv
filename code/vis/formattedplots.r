# O-stats plots with better formatting.
# Author: QDR
# Project: NEON ITV
# Created: 19 Oct 2016
# Last modified: 02 Dec 2016

# Modified 2 Dec: plots with Harvard removed, and work on formatting.
# Modified 7 Nov: Improve scatterplots
# Modified 30 Oct: add continental
# Modified 20 Oct: change axis labels

source('code/vis/loadplotdat.r')

# Jitter plots ------------------------------------------------------------

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

pj <- ggplot(jitterplotdat, aes(x=nullmodel,y=ses)) +
  geom_hline(yintercept=0, linetype='dotted', color = 'blue', size=1) +
  geom_jitter(aes(color=significant), height=0, width=0.25) +
  geom_text(aes(x,y,label=lab), data=jitterplottext, family = 'Helvetica') +
  scale_x_discrete(name =  'Null model', labels = c('Local','Regional')) +
  scale_y_continuous(name = expression(paste('SES of NO'[local]))) +
  scale_color_manual(values = c('gray75', 'black')) +
  theme_john + theme(legend.position = c(0.88,0.1))

ggsave('figs/msfigs/ostat_jitterplot.png', pj, height=5, width=5, dpi=400) 


# Density plots -----------------------------------------------------------



# Local

good_sites <- unique(o2015$siteID)
sites_temporder <- neonsitedata %>% arrange(bio6) %>% dplyr::select(siteID, bio6)

pdensshade <- ggplot(filter(mam_capture_sitemerge, year == 2015, !siteID %in% c('HEAL','DELA','DSNY')) %>% mutate(siteID = factor(siteID, levels=sites_temporder$siteID))) +
  stat_density(adjust = 2, size = 1, aes(x = log10(weight), group = taxonID), fill = 'black', alpha = 0.25, geom='polygon', position = 'identity') + facet_wrap(~ siteID) +
  scale_y_continuous(name = 'probability density', expand = c(0,0)) +
  scale_x_continuous(name = 'body mass (g)', breaks = c(1, 2, 3), labels = c(10, 100, 1000), limits = c(0.5, 3.1)) +
  geom_text(aes(label = paste0('Overlap = ', round(ostat_norm,3)), x = 1.2, y = 13), color = 'black', data = o2015 %>% filter(!siteID %in% c('HEAL','DELA','DSNY'), trait %in% 'logweight') %>% mutate(siteID = factor(siteID, levels=sites_temporder$siteID)), family = 'Helvetica') +
  geom_text(aes(label = paste0('MTCM = ', round(bio6, 1), '°C'), x = 1.2, y = 15), color = 'black', data = neonsitedata %>% filter(siteID %in% good_sites) %>% mutate(siteID = factor(siteID, levels=sites_temporder$siteID)), family = 'Helvetica') +
  theme_john + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), strip.background = element_blank())

source('~/GitHub/NEON/code/facetadjust.r')

png('C:/Users/Q/google_drive/NEON_EAGER/Figures/msfigs2017jan/figs5.png', height = 10, width = 12, res=400, units='in')
  #facetAdjust(pdensshade)
  pdensshade
dev.off()

# Regional

regpoolplotdat <- list()

for (i in 1:length(siteregpoolsp_mam15_iucn)) {
  regpoolplotdat[[i]] <- data.frame(siteID = names(siteregpoollist_mam15_iucn)[i],
                                    taxonID = siteregpoolsp_mam15_iucn[[i]],
                                    logweight = log10(siteregpoollist_mam15_iucn[[i]]$weight))
}
regpoolplotdat <- do.call('rbind',regpoolplotdat)

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

or2015goodsites <- filter(or2015, !siteID %in% c('DELA','DSNY','HEAL'))
or2015goodsites <- or2015goodsites %>%
  mutate(allpool_significant = ostat_reg_allpoolnull_ses < ostat_reg_allpoolnull_seslower | ostat_reg_allpoolnull_ses > ostat_reg_allpoolnull_sesupper,
         bysp_significant = ostat_reg_byspnull_ses < ostat_reg_byspnull_seslower | ostat_reg_byspnull_ses > ostat_reg_byspnull_sesupper)

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

pdenslabels <- ggplot(filter(mam_capture_sitemerge, year == 2015, !siteID %in% c('HEAL','DELA','DSNY')) %>% mutate(siteID = factor(siteID, levels=sites_temporder$siteID), logweight=log10(weight))) +
  stat_density(adjust = 1, size = 1, aes(x = logweight), fill = 'black', alpha = 1, geom = 'polygon', position = 'identity', data = regpoolplotdat %>% filter(!siteID %in% c('HEAL','DELA','DSNY'))) +
  stat_density(adjust = 1, size = 1, aes(x = logweight), fill = 'skyblue', alpha = 0.75, geom='polygon', position = 'identity') + facet_wrap(~ siteID) +
  geom_text(aes(label = paste('NM1:',stattextallpool), x = 2, y = 5), data = or2015regstat %>% filter(!siteID %in% c('DELA','HEAL','DSNY')) %>% mutate(siteID = factor(siteID, levels=sites_temporder$siteID)), family = 'Helvetica') +
  geom_text(aes(label = paste('NM2:',stattextbysp), x = 2, y = 4), data = or2015regstat %>% filter(!siteID %in% c('DELA','HEAL','DSNY')) %>% mutate(siteID = factor(siteID, levels=sites_temporder$siteID)), family = 'Helvetica') +
  geom_text(aes(label = round(ostat_reg,3), x = 2, y = 3), data = or2015regstat %>% filter(!siteID %in% c('DELA','HEAL','DSNY')) %>% mutate(siteID = factor(siteID, levels=sites_temporder$siteID)), family = 'Helvetica') +
  scale_x_continuous(breaks = c(1, 2), labels = c(10, 100), name = expression(paste('log'[10], ' body mass'))) +
  theme_john + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) 
  #qSubtitle('Regional species pools (black) and local communities (blue)', 'Significance of "regional overlap stat" shown, sites ordered by temp')

source('~/GitHub/NEON/code/facetadjust.r')

png('figs/msfigs/regionalspeciespoolsdensity_withlabels.png', height=10, width=12, res=400, units='in')
  facetAdjust(pdenslabels)
dev.off()

pdensconti <-  ggplot(filter(mam_capture_sitemerge, year == 2015, !siteID %in% c('HEAL','DELA','DSNY')) %>% mutate(siteID = factor(siteID, levels=sites_temporder$siteID), logweight=log10(weight))) +
  stat_density(adjust = 1, size = 1, aes(x = logweight), fill = 'black', alpha = 1, geom = 'polygon', position = 'identity', data = filter(mam_capture_sitemerge,  !siteID %in% c('HEAL','DELA','DSNY')) %>% select(-siteID) %>% mutate(logweight=log10(weight))) +
  stat_density(adjust = 1, size = 1, aes(x = logweight), fill = 'skyblue', alpha = 0.75, geom='polygon', position = 'identity') + facet_wrap(~ siteID) +
  geom_text(aes(label = paste('NM1:',stattextallpool), x = 1.5, y = 5), data = orc2015regstat %>% filter(!siteID %in% c('DELA','HEAL','DSNY')) %>% mutate(siteID = factor(siteID, levels=sites_temporder$siteID)), family = 'Helvetica') +
  geom_text(aes(label = paste('NM2:',stattextbysp), x = 1.5, y = 4), data = orc2015regstat %>% filter(!siteID %in% c('DELA','HEAL','DSNY')) %>% mutate(siteID = factor(siteID, levels=sites_temporder$siteID)), family = 'Helvetica') +
  geom_text(aes(label = round(ostat_reg,3), x = 1.5, y = 3), data = orc2015regstat %>% filter(!siteID %in% c('DELA','HEAL','DSNY')) %>% mutate(siteID = factor(siteID, levels=sites_temporder$siteID))) +
  scale_x_continuous(breaks = c(1, 2), labels = c(10, 100), name = expression(paste('log'[10], ' body mass'))) +
  theme_john + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  #qSubtitle('Continental species pools (black) and local communities (blue)', 'Significance of "regional overlap stat" shown, sites ordered by temp')

png('figs/msfigs/continentalspeciespoolsdensity_withlabels.png', height=10, width=12, res=400, units='in')
facetAdjust(pdensconti)
dev.off()


# Scatter plots -----------------------------------------------------------

# Local

porawtemp <- ggplot(o2015 %>% filter(trait=='logweight'), aes(x=bio1)) + 
  geom_segment(aes(y = ostat_norm_localnull_lower, yend = ostat_norm_localnull_upper, xend=bio1), alpha = 0.5, size = 1.5, color = 'skyblue') +
#  geom_segment(aes(y = ostat_norm_regnull_lower, yend = ostat_norm_regnull_upper, xend=bio1), alpha = 0.5, size = 1.5, color = 'goldenrod') +
  geom_point(aes(y = ostat_norm), size = 1.5) +
  labs(y = 'median overlap', x = parse(text=bioclimnames[1])) +
  theme_john
  #qSubtitle('Overlap statistics for 2015 NEON mammals versus MAT', 'Raw values, local and regional nulls')

porawchao <- ggplot(o2015 %>% filter(trait=='logweight'), aes(x=chao1)) + 
  geom_segment(aes(y = ostat_norm_localnull_lower, yend = ostat_norm_localnull_upper, xend=chao1), alpha = 0.5, size = 1.5, color = 'skyblue') +
#  geom_segment(aes(y = ostat_norm_regnull_lower, yend = ostat_norm_regnull_upper, xend=chao1), alpha = 0.5, size = 1.5, color = 'goldenrod') +
  geom_point(aes(y = ostat_norm), size = 1.5) +
  labs(y = 'median overlap', x = 'Species Richness (Chao1)') +
  theme_john
  #qSubtitle('Overlap statistics for 2015 NEON mammals versus Richness', 'Raw values, local and regional nulls')

porawmpd <- ggplot(o2015 %>% filter(trait=='logweight'), aes(x=mpd_z)) + 
  geom_segment(aes(y = ostat_norm_localnull_lower, yend = ostat_norm_localnull_upper, xend=mpd_z), alpha = 0.5, size = 1.5, color = 'skyblue') +
#  geom_segment(aes(y = ostat_norm_regnull_lower, yend = ostat_norm_regnull_upper, xend=mpd_z), alpha = 0.5, size = 1.5, color = 'goldenrod') +
  geom_point(aes(y = ostat_norm), size = 1.5) +
  labs(y = 'median overlap', x = 'Mean Pairwise Distance SES') +
  theme_john
  #qSubtitle('Overlap statistics for 2015 NEON mammals versus MPD', 'Raw values, local and regional nulls')

porawprecip <- ggplot(o2015 %>% filter(trait=='logweight'), aes(x=bio12)) + 
  geom_segment(aes(y = ostat_norm_localnull_lower, yend = ostat_norm_localnull_upper, xend=bio12), alpha = 0.5, size = 1.5, color = 'skyblue') +
#  geom_segment(aes(y = ostat_norm_regnull_lower, yend = ostat_norm_regnull_upper, xend=bio12), alpha = 0.5, size = 1.5, color = 'goldenrod') +
  geom_point(aes(y = ostat_norm), size = 1.5) +
  labs(y = 'median overlap', x = bioclimnames[12]) +
  theme_john

porawtempseas <- ggplot(o2015 %>% filter(trait=='logweight'), aes(x=bio4)) + 
  geom_segment(aes(y = ostat_norm_localnull_lower, yend = ostat_norm_localnull_upper, xend=bio4), alpha = 0.5, size = 1.5, color = 'skyblue') +
#  geom_segment(aes(y = ostat_norm_regnull_lower, yend = ostat_norm_regnull_upper, xend=bio4), alpha = 0.5, size = 1.5, color = 'goldenrod') +
  geom_point(aes(y = ostat_norm), size = 1.5) +
  labs(y = 'median overlap', x = bioclimnames[4]) +
  theme_john

porawprecipseas <- ggplot(o2015 %>% filter(trait=='logweight'), aes(x=bio15)) + 
  geom_segment(aes(y = ostat_norm_localnull_lower, yend = ostat_norm_localnull_upper, xend=bio15), alpha = 0.5, size = 1.5, color = 'skyblue') +
 # geom_segment(aes(y = ostat_norm_regnull_lower, yend = ostat_norm_regnull_upper, xend=bio15), alpha = 0.5, size = 1.5, color = 'goldenrod') +
  geom_point(aes(y = ostat_norm), size = 1.5) +
  labs(y = 'median overlap', x = bioclimnames[15]) +
  theme_john

porawtempcv <- ggplot(o2015 %>% filter(trait=='logweight'), aes(x=cv_bio1)) + 
  geom_segment(aes(y = ostat_norm_localnull_lower, yend = ostat_norm_localnull_upper, xend=cv_bio1), alpha = 0.5, size = 1.5, color = 'skyblue') +
#  geom_segment(aes(y = ostat_norm_regnull_lower, yend = ostat_norm_regnull_upper, xend=cv_bio1), alpha = 0.5, size = 1.5, color = 'goldenrod') +
  geom_point(aes(y = ostat_norm), size = 1.5) +
  labs(y = 'median overlap', x = 'Among-year temperature CV') +
  theme_john

porawprecipcv <- ggplot(o2015 %>% filter(trait=='logweight'), aes(x=cv_bio12)) + 
  geom_segment(aes(y = ostat_norm_localnull_lower, yend = ostat_norm_localnull_upper, xend=cv_bio12), alpha = 0.5, size = 1.5, color = 'skyblue') +
 # geom_segment(aes(y = ostat_norm_regnull_lower, yend = ostat_norm_regnull_upper, xend=cv_bio12), alpha = 0.5, size = 1.5, color = 'goldenrod') +
  geom_point(aes(y = ostat_norm), size = 1.5) +
  labs(y = 'median overlap', x = 'Among-year precipitation CV') +
  theme_john

library(directlabels)
sitelab <- geom_dl(aes(label = siteID, y = ostat_norm), method = list('top.bumptwice', cex = 0.75, vjust = -0.5, fontfamily = 'Helvetica'))
fp <- 'C:/Users/Q/Google Drive/NEON_EAGER/Figures/msfigs2017jan/figs6'
ggsave(file.path(fp,'scatterlocaltemp.png'), porawtemp + sitelab, height=5, width=5, dpi=400)
ggsave(file.path(fp,'scatterlocalprecip.png'), porawprecip + sitelab, height=5, width=5, dpi=400)
ggsave(file.path(fp,'scatterlocalrichness.png'), porawchao + sitelab, height=5, width=5, dpi=400)
ggsave(file.path(fp,'scatterlocalmpd.png'), porawmpd + sitelab, height=5, width=5, dpi=400)
ggsave(file.path(fp,'scatterlocaltempseas.png'), porawtempseas + sitelab, height=5, width=5, dpi=400)
ggsave(file.path(fp,'scatterlocalprecipseas.png'), porawprecipseas + sitelab, height=5, width=5, dpi=400)
ggsave(file.path(fp,'scatterlocaltempcv.png'), porawtempcv + sitelab, height=5, width=5, dpi=400)
ggsave(file.path(fp,'scatterlocalprecipcv.png'), porawprecipcv + sitelab, height=5, width=5, dpi=400)

# Regional

porrawtemp <- ggplot(or2015goodsites %>% filter(trait=='logweight'), aes(x=bio1)) + 
  geom_segment(aes(y = ostat_reg_allpoolnull_lower, yend = ostat_reg_allpoolnull_upper, xend=bio1), alpha = 0.5, size = 1.5, color = 'seagreen3') +
  geom_segment(aes(y = ostat_reg_byspnull_lower, yend = ostat_reg_byspnull_upper, xend=bio1), alpha = 0.5, size = 1.5, color = 'plum2') +
  geom_point(aes(y = ostat_reg), size = 1.5) +
  labs(y = expression(NO[regional]), x = parse(text=bioclimnames[1])) +
  theme_john
  #qSubtitle('Regional O-stats for 2015 NEON mammals versus MAT', 'Raw values, all-pool and by-species nulls')

porrawprecip <- ggplot(or2015goodsites %>% filter(trait=='logweight'), aes(x=bio12)) + 
  geom_segment(aes(y = ostat_reg_allpoolnull_lower, yend = ostat_reg_allpoolnull_upper, xend=bio12), alpha = 0.5, size = 1.5, color = 'seagreen3') +
  geom_segment(aes(y = ostat_reg_byspnull_lower, yend = ostat_reg_byspnull_upper, xend=bio12), alpha = 0.5, size = 1.5, color = 'plum2') +
  geom_point(aes(y = ostat_reg), size = 1.5) +
  labs(y = expression(NO[regional]), x = bioclimnames[12]) +
  theme_john


porrawchao <- ggplot(or2015goodsites %>% filter(trait=='logweight'), aes(x=chao1)) + 
  geom_segment(aes(y = ostat_reg_allpoolnull_lower, yend = ostat_reg_allpoolnull_upper, xend=chao1), alpha = 0.5, size = 1.5, color = 'seagreen3') +
  geom_segment(aes(y = ostat_reg_byspnull_lower, yend = ostat_reg_byspnull_upper, xend=chao1), alpha = 0.5, size = 1.5, color = 'plum2') +
  geom_point(aes(y = ostat_reg), size = 1.5) +
  labs(y = expression(NO[regional]), x = 'Chao1 Species Richness') +
  theme_john
  #qSubtitle('Regional O-stats for 2015 NEON mammals versus richness', 'Raw values, all-pool and by-species nulls')

porrawmpd <- ggplot(or2015goodsites %>% filter(trait=='logweight'), aes(x=mpd_z)) + 
  geom_segment(aes(y = ostat_reg_allpoolnull_lower, yend = ostat_reg_allpoolnull_upper, xend=mpd_z), alpha = 0.5, size = 1.5, color = 'seagreen3') +
  geom_segment(aes(y = ostat_reg_byspnull_lower, yend = ostat_reg_byspnull_upper, xend=mpd_z), alpha = 0.5, size = 1.5, color = 'plum2') +
  geom_point(aes(y = ostat_reg), size = 1.5) +
  labs(y = expression(NO[regional]), x = 'Mean Pairwise Distance SES') +
  theme_john
  #qSubtitle('Regional O-stats for 2015 NEON mammals versus mpd', 'Raw values, all-pool and by-species nulls')

porrawtempseas <- ggplot(or2015goodsites %>% filter(trait=='logweight'), aes(x=bio4)) + 
  geom_segment(aes(y = ostat_reg_allpoolnull_lower, yend = ostat_reg_allpoolnull_upper, xend=bio4), alpha = 0.5, size = 1.5, color = 'seagreen3') +
  geom_segment(aes(y = ostat_reg_byspnull_lower, yend = ostat_reg_byspnull_upper, xend=bio4), alpha = 0.5, size = 1.5, color = 'plum2') +
  geom_point(aes(y = ostat_reg), size = 1.5) +
  labs(y = expression(NO[regional]), x = bioclimnames[4]) +
  theme_john

porrawprecipseas <- ggplot(or2015goodsites %>% filter(trait=='logweight'), aes(x=bio15)) + 
  geom_segment(aes(y = ostat_reg_allpoolnull_lower, yend = ostat_reg_allpoolnull_upper, xend=bio15), alpha = 0.5, size = 1.5, color = 'seagreen3') +
  geom_segment(aes(y = ostat_reg_byspnull_lower, yend = ostat_reg_byspnull_upper, xend=bio15), alpha = 0.5, size = 1.5, color = 'plum2') +
  geom_point(aes(y = ostat_reg), size = 1.5) +
  labs(y = expression(NO[regional]), x = bioclimnames[15]) +
  theme_john

porrawtempcv <- ggplot(or2015goodsites %>% filter(trait=='logweight'), aes(x=cv_bio1)) + 
  geom_segment(aes(y = ostat_reg_allpoolnull_lower, yend = ostat_reg_allpoolnull_upper, xend=cv_bio1), alpha = 0.5, size = 1.5, color = 'seagreen3') +
  geom_segment(aes(y = ostat_reg_byspnull_lower, yend = ostat_reg_byspnull_upper, xend=cv_bio1), alpha = 0.5, size = 1.5, color = 'plum2') +
  geom_point(aes(y = ostat_reg), size = 1.5) +
  labs(y = expression(NO[regional]), x = 'Among-year temperature CV') +
  theme_john

porrawprecipcv <- ggplot(or2015goodsites %>% filter(trait=='logweight'), aes(x=cv_bio12)) + 
  geom_segment(aes(y = ostat_reg_allpoolnull_lower, yend = ostat_reg_allpoolnull_upper, xend=cv_bio12), alpha = 0.5, size = 1.5, color = 'seagreen3') +
  geom_segment(aes(y = ostat_reg_byspnull_lower, yend = ostat_reg_byspnull_upper, xend=cv_bio12), alpha = 0.5, size = 1.5, color = 'plum2') +
  geom_point(aes(y = ostat_reg), size = 1.5) +
  labs(y = expression(NO[regional]), x = 'Among-year precipitation CV') +
  theme_john

sitelabreg <- geom_dl(aes(label = siteID, y = ostat_reg), method = list('top.bumptwice', cex = 0.75, vjust = -0.5, fontfamily = 'Helvetica'))
ggsave('figs/msfigs/scatterregionaltemp.png', porrawtemp + sitelabreg, height=5, width=5, dpi=400)
ggsave('figs/msfigs/scatterregionalprecip.png', porrawprecip + sitelabreg, height=5, width=5, dpi=400)
ggsave('figs/msfigs/scatterregionalrichness.png', porrawchao + sitelabreg, height=5, width=5, dpi=400)
ggsave('figs/msfigs/scatterregionalmpd.png', porrawmpd + sitelabreg, height=5, width=5, dpi=400)
ggsave('figs/msfigs/scatterregionaltempseas.png', porrawtempseas + sitelabreg, height=5, width=5, dpi=400)
ggsave('figs/msfigs/scatterregionalprecipseas.png', porrawprecipseas + sitelabreg, height=5, width=5, dpi=400)
ggsave('figs/msfigs/scatterregionaltempcv.png', porrawtempcv + sitelabreg, height=5, width=5, dpi=400)
ggsave('figs/msfigs/scatterregionalprecipcv.png', porrawprecipcv + sitelabreg, height=5, width=5, dpi=400)


# 27 Oct: simpler plots with logistic line --------------------------------

fx <- function(x, b0, b1) 1/(1 + exp(-(b0 + b1 * x)))
tempco <- summary(reglocalbio)$coeff$mean[,1]
chaoco <- summary(reglocalchao)$coeff$mean[,1]
csc <- scale_color_manual(values = c('gray75','black'))

porawtemp <- ggplot(o2015goodsites %>% filter(trait=='logweight'), aes(x=bio1)) + 
  #geom_segment(aes(y = ostat_norm_localnull_lower, yend = ostat_norm_localnull_upper, xend=bio1), alpha = 0.5, size = 1.5, color = 'skyblue') +
  #geom_segment(aes(y = ostat_norm_regnull_lower, yend = ostat_norm_regnull_upper, xend=bio1), alpha = 0.5, size = 1.5, color = 'goldenrod') +
  geom_point(aes(y = ostat_norm, color=local_significant), size = 2) +
  stat_function(geom='line', fun = fx, args=list(b0 = tempco[1], b1 = tempco[2]), color = 'blue', size = 1.5, n=9999) +
  labs(y = 'Niche Overlap', x = parse(text=bioclimnames[1])) +
  theme_john + theme(legend.position = 'none') + csc
#qSubtitle('Overlap statistics for 2015 NEON mammals versus MAT', 'Raw values, local and regional nulls')

porawchao <- ggplot(o2015goodsites %>% filter(trait=='logweight'), aes(x=chao1)) + 
  #geom_segment(aes(y = ostat_norm_localnull_lower, yend = ostat_norm_localnull_upper, xend=chao1), alpha = 0.5, size = 1.5, color = 'skyblue') +
  #geom_segment(aes(y = ostat_norm_regnull_lower, yend = ostat_norm_regnull_upper, xend=chao1), alpha = 0.5, size = 1.5, color = 'goldenrod') +
  geom_point(aes(y = ostat_norm, color=local_significant), size = 2) +
  stat_function(geom='line', fun = fx, args=list(b0 = chaoco[1], b1 = chaoco[2]), color = 'blue', size = 1.5, n=9999) +
  labs(y = 'Niche Overlap', x = 'Species Richness (Chao1)') +
  theme_john + theme(legend.position = 'none') + csc

porawmpd <- ggplot(o2015goodsites %>% filter(trait=='logweight'), aes(x=mpd_z)) + 
  #geom_segment(aes(y = ostat_norm_localnull_lower, yend = ostat_norm_localnull_upper, xend=mpd_z), alpha = 0.5, size = 1.5, color = 'skyblue') +
  #geom_segment(aes(y = ostat_norm_regnull_lower, yend = ostat_norm_regnull_upper, xend=mpd_z), alpha = 0.5, size = 1.5, color = 'goldenrod') +
  geom_point(aes(y = ostat_norm, color=local_significant), size = 2) +
  labs(y = 'Niche Overlap', x = 'Mean Pairwise Distance SES') +
  theme_john+ theme(legend.position = 'none') + csc

ggsave('figs/msfigs/simplescatterlocaltemp.png', porawtemp + theme(aspect.ratio=1), height=5, width=5, dpi=400)
ggsave('figs/msfigs/simplescatterlocalrichness.png', porawchao + theme(aspect.ratio=1), height=5, width=5, dpi=400)

# Make these into figures a and b

panela <- porawtemp# + theme(aspect.ratio=1)
panelb <- porawchao + theme(axis.text.y=element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank())
panelc <- porawmpd + theme(axis.text.y=element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank())


library(gridExtra)
grid.arrange(panela, panelb, nrow=1, widths=c(1.05, 0.95))

library(grid)
png('figs/msfigs/simplescatter3panels.png', height=4, width=12, res=400, units='in')
grid.newpage()
grid.draw(cbind(ggplotGrob(panela), ggplotGrob(panelb), ggplotGrob(panelc), size = "last"))
dev.off()


# Logit scale scatterplots ------------------------------------------------

library(scales)
inverse_logit_trans <- trans_new("inverse logit",
                                 transform = plogis,
                                 inverse = qlogis)

logit_trans <- trans_new("logit",
                                 transform = qlogis,
                                 inverse = plogis)

porawtemp + scale_y_continuous(trans = inverse_logit_trans)

porawtemp <- ggplot(o2015goodsites %>% filter(trait=='logweight'), aes(x=bio1)) + 
  geom_point(aes(y = ostat_norm, color=local_significant), size = 3) +
 # geom_abline(intercept = tempco[1], slope = tempco[2], color = 'dodgerblue', size = 1.5) +
  stat_function(geom='line', fun = fx, args=list(b0 = tempco[1], b1 = tempco[2]), color = 'dodgerblue', size = 1.5, n=9999) +
  labs(y = 'Niche Overlap', x = parse(text=bioclimnames[1])) +
  theme_john + theme(legend.position = 'none') + csc +
  coord_trans(y = logit_trans) +
  scale_x_continuous(expand = c(0,0), breaks = c(0,10,20), labels=c(0,10,20), limits=c(-0.5,21.5))

porawchao <- ggplot(o2015goodsites %>% filter(trait=='logweight'), aes(x=chao1)) + 
  geom_point(aes(y = ostat_norm, color=local_significant), size = 3) +
  #geom_abline(intercept = chaoco[1], slope = chaoco[2], color = 'dodgerblue', size = 1.5) +
  stat_function(geom='line', fun = fx, args=list(b0 = chaoco[1], b1 = chaoco[2]), color = 'dodgerblue', size = 1.5, n=9999) +
  labs(y = 'Niche Overlap', x = 'Species Richness (Chao1)') +
  theme_john + theme(legend.position = 'none') + csc +
  coord_trans(y = logit_trans) +
  scale_x_continuous(expand = c(0,0), breaks = c(5,10,15), labels=c(5,10,15), limits=c(4,16))

# Regional

mpdco <- summary(regregionalmpd)$coeff$mean[,1]
chaoco <- summary(regregionalchao)$coeff$mean[,1]

porrawtemp <- ggplot(or2015goodsites %>% filter(trait=='logweight'), aes(x=bio1)) + 
 # geom_segment(aes(y = ostat_reg_allpoolnull_lower, yend = ostat_reg_allpoolnull_upper, xend=bio1), alpha = 0.5, size = 1.5, color = 'seagreen3') +
#  geom_segment(aes(y = ostat_reg_byspnull_lower, yend = ostat_reg_byspnull_upper, xend=bio1), alpha = 0.5, size = 1.5, color = 'plum2') +
  geom_point(aes(y = ostat_reg), size = 1.5) +
  labs(y = expression(NO[regional]), x = parse(text=bioclimnames[1])) +
  theme_john



porrawchao <- ggplot(or2015goodsites %>% filter(trait=='logweight'), aes(x=chao1)) + 
 # geom_segment(aes(y = ostat_reg_allpoolnull_lower, yend = ostat_reg_allpoolnull_upper, xend=chao1), alpha = 0.5, size = 1.5, color = 'seagreen3') +
#  geom_segment(aes(y = ostat_reg_byspnull_lower, yend = ostat_reg_byspnull_upper, xend=chao1), alpha = 0.5, size = 1.5, color = 'plum2') +
  geom_point(aes(y = ostat_reg), size = 1.5) +
  stat_function(geom='line', fun = fx, args=list(b0 = chaoco[1], b1 = chaoco[2]), color = 'blue', size = 1.5) +
  labs(y = expression(NO[regional]), x = 'Chao1 Species Richness') +
  theme_john
#qSubtitle('Regional O-stats for 2015 NEON mammals versus richness', 'Raw values, all-pool and by-species nulls')

porrawmpd <- ggplot(or2015goodsites %>% filter(trait=='logweight'), aes(x=mpd_z)) + 
 # geom_segment(aes(y = ostat_reg_allpoolnull_lower, yend = ostat_reg_allpoolnull_upper, xend=mpd_z), alpha = 0.5, size = 1.5, color = 'seagreen3') +
#  geom_segment(aes(y = ostat_reg_byspnull_lower, yend = ostat_reg_byspnull_upper, xend=mpd_z), alpha = 0.5, size = 1.5, color = 'plum2') +
  geom_point(aes(y = ostat_reg), size = 1.5) +
  stat_function(geom='line', fun = fx, args=list(b0 = mpdco[1], b1 = mpdco[2]), color = 'blue', size = 1.5) +
  labs(y = expression(NO[regional]), x = 'Mean Pairwise Distance SES') +
  theme_john

panela <- porrawtemp# + theme(aspect.ratio=1)
panelb <- porrawchao + theme(axis.text.y=element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank())
panelc <- porrawmpd + theme(axis.text.y=element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank())


library(grid)
png('figs/msfigs/simplescatter3panels_regional.png', height=5, width=15, res=400, units='in')
grid.newpage()
grid.draw(cbind(ggplotGrob(panela), ggplotGrob(panelb), ggplotGrob(panelc), size = "last"))
dev.off()


# Continental
pocrawtemp <- ggplot(orc2015goodsites %>% filter(trait=='logweight'), aes(x=bio1)) + 
  # geom_segment(aes(y = ostat_reg_allpoolnull_lower, yend = ostat_reg_allpoolnull_upper, xend=bio1), alpha = 0.5, size = 1.5, color = 'seagreen3') +
  #  geom_segment(aes(y = ostat_reg_byspnull_lower, yend = ostat_reg_byspnull_upper, xend=bio1), alpha = 0.5, size = 1.5, color = 'plum2') +
  geom_point(aes(y = ostat_reg), size = 1.5) +
  labs(y = expression(NO[regional]), x = parse(text=bioclimnames[1])) +
  theme_john

pocrawchao <- ggplot(orc2015goodsites %>% filter(trait=='logweight'), aes(x=chao1)) + 
  # geom_segment(aes(y = ostat_reg_allpoolnull_lower, yend = ostat_reg_allpoolnull_upper, xend=chao1), alpha = 0.5, size = 1.5, color = 'seagreen3') +
  #  geom_segment(aes(y = ostat_reg_byspnull_lower, yend = ostat_reg_byspnull_upper, xend=chao1), alpha = 0.5, size = 1.5, color = 'plum2') +
  geom_point(aes(y = ostat_reg), size = 1.5) +
  #stat_function(geom='line', fun = fx, args=list(b0 = chaoco[1], b1 = chaoco[2]), color = 'blue', size = 1.5) +
  labs(y = expression(NO[regional]), x = 'Chao1 Species Richness') +
  theme_john
#qSubtitle('Regional O-stats for 2015 NEON mammals versus richness', 'Raw values, all-pool and by-species nulls')

pocrawmpd <- ggplot(orc2015goodsites %>% filter(trait=='logweight'), aes(x=mpd_z)) + 
  # geom_segment(aes(y = ostat_reg_allpoolnull_lower, yend = ostat_reg_allpoolnull_upper, xend=mpd_z), alpha = 0.5, size = 1.5, color = 'seagreen3') +
  #  geom_segment(aes(y = ostat_reg_byspnull_lower, yend = ostat_reg_byspnull_upper, xend=mpd_z), alpha = 0.5, size = 1.5, color = 'plum2') +
  geom_point(aes(y = ostat_reg), size = 1.5) +
 # stat_function(geom='line', fun = fx, args=list(b0 = mpdco[1], b1 = mpdco[2]), color = 'blue', size = 1.5) +
  labs(y = expression(NO[regional]), x = 'Mean Pairwise Distance SES') +
  theme_john

panela <- pocrawtemp# + theme(aspect.ratio=1)
panelb <- pocrawchao + theme(axis.text.y=element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank())
panelc <- pocrawmpd + theme(axis.text.y=element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank())


library(grid)
png('figs/msfigs/simplescatter3panels_continental.png', height=5, width=15, res=400, units='in')
grid.newpage()
grid.draw(cbind(ggplotGrob(panela), ggplotGrob(panelb), ggplotGrob(panelc), size = "last"))
dev.off()
