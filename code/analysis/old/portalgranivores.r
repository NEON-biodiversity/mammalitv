# Portal Rodent Data
# See http://esapubs.org/archive/ecol/E090/118/metadata.htm
# Last modified: 19 Oct 2016

# Modified 20 Apr 2018: specify data paths to hpcc

data_path <- '/mnt/research/neon'

# Added 19 Oct: each pairwise overlap from each year
# Added 17 Oct: climate time series
# Added 11 Oct: time series analysis

library(dplyr)
library(ggplot2)

portal <- read.csv(file.path(data_path, 'external_data/raw_external_data/other/Portal_rodents_19772002.csv'))

spcodes <- c('DM','DO','DS','PF')

portal4 <- portal %>% filter(species %in% spcodes, !is.na(wgt))

# View density plots for the four species

source('code/analysis/densityoverlap.r')

portaloverlap4 <- community_overlap_wmedian(traits = log10(portal4$wgt), sp = as.character(portal4$species))
portaloverlapall <- community_overlap_wmedian(traits = log10(portal$wgt), sp = as.character(portal$species))

# Null distr. for portal granivores
portaloverlap_null <- numeric(999)
for (i in 1:999) portaloverlap_null[i] <- community_overlap_wmedian(traits = log10(portal4$wgt), sp = sample(as.character(portal4$species)))
quantile(portaloverlap_null, probs=c(0.025,0.975))

portaloverlapall_null <- numeric(999)
for (i in 1:999) portaloverlapall_null[i] <- community_overlap_wmedian(traits = log10(portal$wgt), sp = sample(as.character(portal$species)))
quantile(portaloverlapall_null, probs=c(0.025,0.975))

ggplot(portal4, aes(x = log10(wgt), group = species)) + geom_density(aes(fill=species), alpha = 0.25) + theme_john +
  geom_text(aes(label = round(portaloverlap4, 3)), x = Inf, y = Inf, hjust = 1, vjust = 1) +
  ggtitle('Portal Granivores')
ggsave('figs/png/tstats/portaloverlapgranivores.png', height=5, width=6, dpi=400)

ggplot(portal, aes(x = log10(wgt), group = species)) + geom_density(aes(fill=species), alpha = 0.25) + theme_john +
  geom_text(aes(label = round(portaloverlapall, 3)), x = Inf, y = Inf, hjust = 1, vjust = 1) +
  ggtitle('Portal All Species')
ggsave('figs/png/tstats/portaloverlapallspecies.png', height=5, width=7, dpi=400)


# Portal overlap over time (added 11 Oct) ---------------------------------

source('code/analysis/densityoverlap.r')
library(dplyr)
portaloverlap4_time <- portal4 %>%
  group_by(yr) %>%
  summarize(overlap = community_overlap_wmedian(traits = log10(wgt), sp = as.character(species)),
            n = n(),
            S = length(unique(species)))

plot(overlap~yr, data=portaloverlap4_time)
pairs(portaloverlap4_time)

# Auto-regressive model
portalar1 <- ar(x = ts(portaloverlap4_time$overlap), aic=T, order.max=1)

# Highest overlap year was 2000.
ggplot(portal4 %>% filter(yr==2000), aes(x = log10(wgt), group = species)) + geom_density(aes(fill=species), alpha = 0.25) + theme_bw() +
  #geom_text(aes(label = round(portaloverlap4, 3)), x = Inf, y = Inf, hjust = 1, vjust = 1) +
  ggtitle('Portal Granivores 2000')

ggplot(portal4 %>% filter(yr==2001), aes(x = log10(wgt), group = species)) + geom_density(aes(fill=species), alpha = 0.25) + theme_bw() +
  #geom_text(aes(label = round(portaloverlap4, 3)), x = Inf, y = Inf, hjust = 1, vjust = 1) +
  ggtitle('Portal Granivores 2001')

ggplot(portal4 %>% filter(yr==2002), aes(x = log10(wgt), group = species)) + geom_density(aes(fill=species), alpha = 0.25) + theme_bw() +
  #geom_text(aes(label = round(portaloverlap4, 3)), x = Inf, y = Inf, hjust = 1, vjust = 1) +
  ggtitle('Portal Granivores 2002')

ggplot(portal4 %>% filter(yr==1990), aes(x = log10(wgt), group = species)) + geom_density(aes(fill=species), alpha = 0.25) + theme_bw() +
  #geom_text(aes(label = round(portaloverlap4, 3)), x = Inf, y = Inf, hjust = 1, vjust = 1) +
  ggtitle('Portal Granivores 1990')


# Climate over time -------------------------------------------------------

load(file.path(data_path, 'MS1_RodentOverlap/R_data/portalprismtempobj.r'))

# Load the climate data into a single data frame.
extr2df <- function(x, name) {
  dates <- sapply(strsplit(dimnames(x)[[2]], '_'), '[', 5)
  months <- as.numeric(substr(dates,5,6))
  years <- as.numeric(substr(dates,1,4))
  res <- data.frame(year = years, month = months, y = as.numeric(x[1,]))
  names(res)[3] <- name
  res
}

pptdat <- extr2df(ppt_extract, 'ppt')
tmeandat <- extr2df(tmean_extract, 'temp')
portalclimate <- full_join(tmeandat, pptdat)
portalclimatebyyear <- portalclimate %>% group_by(year) %>% summarize(ppt = sum(ppt), temp = mean(temp))

portaloverlap4_time <- left_join(portaloverlap4_time, portalclimatebyyear %>% rename(yr = year))


# Model with time and covariates
library(gam)
timetrend <- gam(overlap ~ s(yr), data = portaloverlap4_time)
fittores <- lm(timetrend$res ~ ppt + temp + n, data = portaloverlap4_time)
# No predictive power after fitting time trend
fitraw <- lm(overlap ~ ppt + temp + n, data = portaloverlap4_time) # Possibly a time effect with higher overlap in higher temp.

# Plot all years
portaldensyr <- ggplot(portal4, aes(x = log10(wgt), group = species)) + geom_density(aes(fill=species), alpha = 0.25) + theme_bw() +
  facet_wrap(~ yr) +
  #geom_text(aes(label = round(portaloverlap4, 3)), x = Inf, y = Inf, hjust = 1, vjust = 1) +
  ggtitle('Portal Granivores By Year') +
  theme(panel.grid=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank(), legend.position = 'bottom') +
  scale_x_continuous(breaks = c(1,2), labels=c(10,100), name=expression(paste('log'[10], ' body mass')))

source('code/facetadjust.r')

pdf('figs/pdf/portaloverlapplots.pdf', height=10, width=12)
  facetAdjust(portaldensyr)
dev.off()


# Individual pairwise overlaps for each year ------------------------------

community_overlap_eachpair <- function(traits, sp, norm = TRUE, bw = NULL, n = NULL) {
  sp <- as.character(sp)
  dat <- data.frame(traits=traits, sp=sp, stringsAsFactors = FALSE)
  dat <- dat[complete.cases(dat), ]
  abunds <- table(dat$sp)
  dat <- dat[dat$sp %in% names(abunds)[abunds>1], ]
  traitlist <- split(dat$traits, dat$sp)
  nspp <- length(traitlist)
  
  if (nspp < 2) return(NA)
  
  overlaps <- numeric(0)
  spp1 <- character(0)
  spp2 <- character(0)
  
  for (sp_a in 1:(nspp-1)) {
    for (sp_b in (sp_a+1):nspp) {
      o <- pairwise_overlap(a = traitlist[[sp_a]], b = traitlist[[sp_b]], norm = norm, bw = bw, n = n)
      overlaps <- c(overlaps, o[1])
      spp1 <- c(spp1, names(traitlist)[sp_a])
      spp2 <- c(spp2, names(traitlist)[sp_b])
    }
  }
  
  data.frame(sp_a = spp1, sp_b = spp2, overlap = overlaps)
  
}


portaloverlap4_allpairs <- portal4 %>%
  group_by(yr) %>%
  do(community_overlap_eachpair(traits = log10(.$wgt), sp = as.character(.$species)))

theme_john <- theme_bw() + theme(panel.grid = element_blank(), 
                                 axis.text = element_text(size = 12), 
                                 axis.title = element_text(size = 18),
                                 text = element_text(family = 'Helvetica'))


pcolorlines <- ggplot(portaloverlap4_allpairs %>% mutate(pair = factor(paste(sp_a, sp_b, sep = ':'))), aes(x = yr, y = overlap)) +
  geom_point(aes(color = pair, group = pair)) +
  geom_line(data = portaloverlap4_time, color = 'black', size = 1.5) +
  geom_text(data = portaloverlap4_time, aes(label = S, y = overlap + 0.02), family = 'Helvetica') +
  labs(x = 'Year', color = 'Species pair') +
  theme_john

ppanellines <- ggplot(portaloverlap4_allpairs %>% mutate(pair = factor(paste(sp_a, sp_b, sep = ':'))), aes(x = yr, y = overlap)) +
  geom_line() + facet_wrap(~ pair, scales = 'free_y') + theme_john + labs(x = 'Year')

ggsave('figs/png/portalpairwise1.png', pcolorlines, height=5, width=6, dpi=400)
ggsave('figs/png/portalpairwise2.png', ppanellines, height = 7, width=10, dpi=400)

# Run just the granivores for the neon sites
granoverlapbysite <- mam_capture_sitemerge %>% 
  filter(year==2015, siteID != 'DSNY', Pineda_Main_food=='Granivore') %>% # Expunge Disney because of poor sampling there 
  mutate(logweight = log10(weight)) %>% 
  group_by(siteID) %>%
  do(overlap_norm = community_overlap(traits = .$logweight, sp = .$taxonID, norm = TRUE),
     overlap_unnorm = community_overlap(traits = .$logweight, sp = .$taxonID, norm = FALSE)) 

granoverlapbysite <- with(granoverlapbysite, data.frame(siteID=siteID, overlap_norm=unlist(overlap_norm), overlap_unnorm=unlist(overlap_unnorm)))

ggplot(filter(mam_capture_sitemerge, year == 2015, !siteID %in% c('HEAL','DELA','DSNY'), Pineda_Main_food=='Granivore') %>% mutate(siteID = factor(siteID, levels=sites_temporder$siteID))) +
  stat_density(adjust = 2, size = 1, aes(x = log10(weight), group = taxonID), fill = 'black', alpha = 0.25, geom='polygon', position = 'identity') + facet_wrap(~ siteID) +
  scale_x_continuous(breaks = c(1, 2), labels = c(10, 100)) +
 # geom_text(aes(label = round(overlap_norm,3), x = 2.5, y = 15), color = 'black', data = overlapstatsbysite %>% filter(!siteID %in% c('HEAL','DELA','DSNY')) %>% mutate(siteID = factor(siteID, levels=sites_temporder$siteID))) +
#  geom_text(aes(label = paste(round(bio1, 2), 'C'), x = 1, y = 15), color = 'black', data = neonsitedata %>% filter(siteID %in% good_sites) %>% mutate(siteID = factor(siteID, levels=sites_temporder$siteID))) +
  theme_john + ggtitle('Mammal Weights ordered by temperature')