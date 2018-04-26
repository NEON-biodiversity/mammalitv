# Run to load all data for plotting.
# Created: 19 Oct 2016
# Modified: 27 Oct 2016 (add heterogeneity)
# Modified 02 Dec 2016 (remove Harvard)
# Modified 15 Dec 2016 (add local pool)

rm(list=ls())

# Weighted medians: load files
files2load <- dir('C:/Users/Q/Dropbox/neon/code/overlapstat', pattern='*wmed.r')
for (i in file.path('C:/Users/Q/Dropbox/neon/code/overlapstat', files2load)) load(i)

# Mean-only analysis: load files
load('C:/Users/Q/Dropbox/neon/code/overlapstat/Ostats_bysite2015_noitv.r')

load('C:/Users/Q/Dropbox/neon/code/overlapstat/Ostats_allregional.r')

load('C:/Users/Q/Dropbox/neon/code/overlapstat/Ostats_bysite2015_wmedlocal.r')

#load('C:/Users/Q/Dropbox/neon/code/overlapstat/Ostats_continentpool.r')

source('code/analysis/ostat2longform.r')

o2014 <- ostat2longform(Ostats_bysite2014)
o2015 <- ostat2longform(Ostats_bysite2015)
oallyears <- ostat2longform(Ostats_bysiteallyears)

or2014 <- regionalostat2longform(Ostats_regional2014)
or2015 <- regionalostat2longform(Ostats_regional2015)
orallyears <- regionalostat2longform(Ostats_regionalallyears)

orc2015 <- regionalostat2longform(Ostats_regcontinent2015)

orlocal2015 <- ostat2longform(Ostats_bysite2015local)

noitv2015 <- ostat2longform(Ostats_bysite2015_noITV)

library(dplyr)
library(ggplot2)
library(lubridate)
library(reshape2)
library(extrafont)

source('~/qutil.r')
source('code/bioclimnames.r')

load('C:/Users/Q/Dropbox/neon/code/mammalPDbyplotobject.r')

source('code/data_extraction/loadmammallocal.r')
source('code/data_extraction/createregpoollists.r')

spatialstat <- read.csv('C:/Users/Q/Dropbox/neon/data/external_datasets/NEON_spatial_stats_highres.csv')
spstatmeans <- neonplotdata %>% cbind(spatialstat) %>% group_by(siteID) %>% summarize(ruggedness=mean(tri, na.rm=TRUE))

richnessdf <- read.csv('richness.csv', stringsAsFactors = FALSE)
heterodf <- read.csv('C:/Users/Q/Dropbox/neon/data/heterogeneity.csv', stringsAsFactors = FALSE)

mammalPDsite <- mammalPD %>% mutate(siteID = substr(plotID,1,4)) %>% group_by(siteID) %>% summarize(mpd_z = mean(mpd.obs.z2015, na.rm=T), mntd_z = mean(mntd.obs.z2015, na.rm=T))

o2014 <- o2014 %>% rename(siteID=site) %>% left_join(neonsitedata) %>% left_join(spstatmeans) %>% left_join(mammalPDsite) %>% left_join(richnessdf) %>% left_join(heterodf) %>% filter(chao1 > 1)
o2015 <- o2015 %>% rename(siteID=site) %>% left_join(neonsitedata) %>% left_join(spstatmeans) %>% left_join(mammalPDsite) %>% left_join(richnessdf) %>% left_join(heterodf) %>% filter(chao1 > 1)
oallyears <- oallyears %>% rename(siteID=site) %>% left_join(neonsitedata) %>% left_join(spstatmeans) %>% left_join(mammalPDsite) %>% left_join(richnessdf) %>% left_join(heterodf) %>% filter(chao1 > 1)

or2014 <- or2014 %>% rename(siteID=site) %>% left_join(neonsitedata) %>% left_join(spstatmeans) %>% left_join(mammalPDsite) %>% left_join(richnessdf) %>% left_join(heterodf) %>% filter(chao1 > 1)
or2015 <- or2015 %>% rename(siteID=site) %>% left_join(neonsitedata) %>% left_join(spstatmeans) %>% left_join(mammalPDsite) %>% left_join(richnessdf) %>% left_join(heterodf) %>% filter(chao1 > 1)
orallyears <- orallyears %>% rename(siteID=site) %>% left_join(neonsitedata) %>% left_join(spstatmeans) %>% left_join(mammalPDsite) %>% left_join(richnessdf) %>% left_join(heterodf) %>% filter(chao1 > 1)

orc2015 <- orc2015 %>% rename(siteID=site) %>% left_join(neonsitedata) %>% left_join(spstatmeans) %>% left_join(mammalPDsite) %>% left_join(richnessdf) %>% left_join(heterodf) %>% filter(chao1 > 1)
#orcallyears <- orcallyears %>% rename(siteID=site) %>% left_join(neonsitedata) %>% left_join(spstatmeans) %>% left_join(mammalPDsite) %>% left_join(richnessdf) %>% filter(chao1 > 1)

orlocal2015 <- orlocal2015 %>% rename(siteID=site) %>% left_join(neonsitedata) %>% left_join(spstatmeans) %>% left_join(mammalPDsite) %>% left_join(richnessdf) %>% left_join(heterodf) %>% filter(chao1 > 1)

noitv2015 <- noitv2015 %>% rename(siteID=site) %>% left_join(neonsitedata) %>% left_join(spstatmeans) %>% left_join(mammalPDsite) %>% left_join(richnessdf) %>% left_join(heterodf) %>% filter(chao1 > 1)

o2015goodsites <- filter(o2015, !siteID %in% c('DELA','DSNY','HEAL', 'HARV'))
noitv2015goodsites <- filter(noitv2015, !siteID %in% c('DELA','DSNY','HEAL', 'HARV'))

# Set plot parameters
theme_john <- theme_bw() + theme(panel.grid = element_blank(), 
                                 axis.text = element_text(size = 12), 
                                 axis.title = element_text(size = 18),
                                 text = element_text(family = 'Helvetica'))

# Add some additional parameters to get rid of the axis ticks and numbers if you want.
theme_noaxisnumbers <- theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
