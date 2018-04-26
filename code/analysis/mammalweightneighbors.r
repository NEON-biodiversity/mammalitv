# Mammal trait variance and neighbor distance analysis/visualization
# Author: QDR
# Project: NEON ITV
# Created: 22 Sep 2016
# Last modified: 20 Apr 2018

# Modified 20 Apr 2018: specify paths

# Load and clean data -----------------------------------------------------


# Use the same data cleaning as in the tstat_by_guild.r script to start out with.

data_path <- '/mnt/research/neon'

iucn <- read.csv(file.path(data_path, 'external_data/final_external_data/IUCN_mammal_ranges.csv'))
load(file.path(data_path, 'final_data/allorganismal_latest.r'))
load(file.path(data_path, 'final_data/allsiteplot_latest.r'))

mammalTax <- read.csv(file.path(data_path,'raw_data/NEON_mam_taxonomy.csv'))
mammalTraits <- read.csv(file.path(data_path,'external_data/final_external_data/NEON_miscmammaltraits.csv'))

# Clean mammal data to get just adults.
# Take the median value of the traits for each individual if it was recaptured
library(dplyr)
library(lubridate)
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

mam_noID <- mam_capture[mam_capture$individualandtag=='', c('year', 'siteID','plotID', 'taxonID','family', 'individualandtag','hindfootLength','earLength','tailLength','totalLength','weight','sex','lifeStage','Pineda_Main_food')]

mam_grp <- filter(mam_capture, individualandtag != '') %>%
  group_by(year, siteID, plotID, taxonID, family, individualandtag)

mam_byindiv <- summarize_at(mam_grp, vars(hindfootLength, earLength, tailLength, totalLength, weight), median, na.rm=TRUE)
mam_byindiv_class <- do(mam_grp, data.frame(sex = .$sex[1], lifeStage = names(sort(table(.$lifeStage),decreasing=TRUE))[1], Pineda_Main_food=.$Pineda_Main_food[1]))

mam_byindiv <- cbind(as.data.frame(mam_byindiv), as.data.frame(mam_byindiv_class[,c('sex', 'lifeStage','Pineda_Main_food')]))                  
mam_byindiv <-rbind(mam_byindiv, mam_noID)

mam_capture_sitemerge <- left_join(mam_byindiv, filter(neonplotdata, subtype == 'mammalGrid'))


# Density plot visualizations ---------------------------------------------

library(ggplot2)

# Density plots of mammal sizes (rodents only, no shrews)

ggplot(filter(mam_capture_sitemerge, year == 2015), aes(x = log10(weight), group = taxonID, color = Pineda_Main_food)) +
  geom_density() + facet_wrap(~ plotID)

ggsave('figs/pdf/alldistributions.pdf', height = 16, width = 16)

ggplot(filter(mam_capture_sitemerge, year == 2015), aes(x = log10(weight), group = taxonID, color = Pineda_Main_food)) +
  geom_density() + facet_wrap(~ plotID)


# Pool by site since there are way too many density plots to visualize any pattern with subplot.
ggplot(filter(mam_capture_sitemerge, year == 2015), aes(x = log10(weight), group = taxonID, color = Pineda_Main_food)) +
  geom_density(adjust = 2) + facet_wrap(~ siteID) +
  theme_bw()

# Sort sites by richness of each guild to order the density plots.

# Get richness of each guild at each plot
guildrichness <- mam_capture_sitemerge %>% 
  group_by(siteID, Pineda_Main_food) %>%
  summarize(richness = length(unique(taxonID)))

mam_capture_sitemerge <- left_join(mam_capture_sitemerge, guildrichness)  

plots_granivoreorder <- guildrichness %>% filter(Pineda_Main_food == 'Granivore') %>% arrange(richness)
plots_ggorder <- guildrichness %>% filter(Pineda_Main_food %in% c('Generalist','Granivore')) %>% group_by(siteID) %>% summarize(richness=sum(richness)) %>% arrange(richness)

sites_temporder <- neonsitedata %>% arrange(bio1) %>% dplyr::select(siteID, bio1)

theme_john <- theme_bw() + theme(panel.grid = element_blank(), axis.text = element_text(size = 12), axis.title = element_text(size = 18))

# Granivores only
ggplot(filter(mam_capture_sitemerge, year == 2015, Pineda_Main_food == 'Granivore', richness>1) %>% mutate(siteID = factor(siteID, levels=plots_granivoreorder$siteID)), aes(x = log10(weight), group = taxonID, color = taxonID)) +
  geom_density(adjust = 2, size = 1) + facet_wrap(~ siteID) +
  theme_john + ggtitle('Granivores')

# Generalists + Granivores
# As polygon
ggplot(filter(mam_capture_sitemerge, year == 2015, Pineda_Main_food %in% c('Generalist','Granivore'), richness>1) %>% mutate(siteID = factor(siteID, levels=plots_ggorder$siteID))) +
  geom_density(adjust = 2, size = 1, aes(x = log10(weight), group = taxonID, color = taxonID)) + facet_wrap(~ siteID) +
  geom_text(aes(label = richness, x = 2.5, y = 9), color = 'black', data = plots_ggorder %>% mutate(siteID = factor(siteID, levels=plots_ggorder$siteID))) + 
  theme_john + ggtitle('Generalists + Granivores')
ggsave('figs/png/tstats/densityplotspolygons.png', height = 12, width = 14)

# As line
ggplot(filter(mam_capture_sitemerge, year == 2015, Pineda_Main_food %in% c('Generalist','Granivore'), richness>1) %>% mutate(siteID = factor(siteID, levels=plots_ggorder$siteID))) +
  stat_density(adjust = 2, size = 1, aes(x = log10(weight), group = taxonID, color = taxonID), geom='line', position = 'identity') + facet_wrap(~ siteID) +
  geom_text(aes(label = richness, x = 2.5, y = 9), color = 'black', data = plots_ggorder %>% mutate(siteID = factor(siteID, levels=plots_ggorder$siteID))) + 
  geom_text(aes(label = t, x = 1, y = 9), color = 'black', data = tmgg %>% filter(stat == 'T_IP.IC', trait == 'logweight') %>% group_by(siteID) %>% summarize(t = round(mean(t, na.rm=T), 3)) %>% mutate(siteID = factor(siteID, levels=plots_ggorder$siteID))) +
  theme_john + ggtitle('Generalists + Granivores')
ggsave('figs/png/tstats/densityplotslines.png', height = 12, width = 14)

# As line, sorted by temperature
# Also include overlap stats

good_sites <- overlapstatsbysite$siteID[!is.na(overlapstatsbysite$overlap_norm)]

ggplot(filter(mam_capture_sitemerge, year == 2015, !siteID %in% c('HEAL','DELA','DSNY')) %>% mutate(siteID = factor(siteID, levels=sites_temporder$siteID))) +
  stat_density(adjust = 2, size = 1, aes(x = log10(weight), group = taxonID, color = taxonID), geom='line', position = 'identity') + facet_wrap(~ siteID) +
  scale_x_continuous(breaks = c(1, 2), labels = c(10, 100)) +
  geom_text(aes(label = round(overlap_norm,3), x = 2.5, y = 9), color = 'black', data = overlapstatsbysite %>% filter(!siteID %in% c('HEAL','DELA','DSNY')) %>% mutate(siteID = factor(siteID, levels=sites_temporder$siteID))) +
  geom_text(aes(label = paste(round(bio1, 2), 'C'), x = 1, y = 9), color = 'black', data = neonsitedata %>% filter(siteID %in% good_sites) %>% mutate(siteID = factor(siteID, levels=sites_temporder$siteID))) +
  theme_john + ggtitle('Mammal Weights ordered by temperature')
ggsave('figs/png/tstats/densityplotsbytemperature_lines.png', height = 12, width = 14)

# As shaded partially transparent polygons
ggplot(filter(mam_capture_sitemerge, year == 2015, !siteID %in% c('HEAL','DELA','DSNY')) %>% mutate(siteID = factor(siteID, levels=sites_temporder$siteID))) +
  stat_density(adjust = 2, size = 1, aes(x = log10(weight), group = taxonID), fill = 'black', alpha = 0.25, geom='polygon', position = 'identity') + facet_wrap(~ siteID) +
  scale_x_continuous(breaks = c(1, 2), labels = c(10, 100)) +
  geom_text(aes(label = round(overlap_norm,3), x = 2.5, y = 15), color = 'black', data = overlapstatsbysite %>% filter(!siteID %in% c('HEAL','DELA','DSNY')) %>% mutate(siteID = factor(siteID, levels=sites_temporder$siteID))) +
  geom_text(aes(label = paste(round(bio1, 2), 'C'), x = 1, y = 15), color = 'black', data = neonsitedata %>% filter(siteID %in% good_sites) %>% mutate(siteID = factor(siteID, levels=sites_temporder$siteID))) +
  theme_john + ggtitle('Mammal Weights ordered by temperature')
ggsave('figs/png/tstats/densityplotsbytemperature_shaded.png', height = 10, width = 12)

# As histograms instead of density plots (non normalized)
ggplot(filter(mam_capture_sitemerge, year == 2015, !siteID %in% c('HEAL','DELA','DSNY')) %>% mutate(siteID = factor(siteID, levels=sites_temporder$siteID))) +
  geom_histogram(aes(x = log10(weight), group = taxonID), fill = 'black', alpha = 0.25, position = 'identity') + facet_wrap(~ siteID) +
  scale_x_continuous(breaks = c(1, 2), labels = c(10, 100)) +
  geom_text(aes(label = round(overlap_norm,3), x = 2.5, y = 15), color = 'black', data = overlapstatsbysite %>% filter(!siteID %in% c('HEAL','DELA','DSNY')) %>% mutate(siteID = factor(siteID, levels=sites_temporder$siteID))) +
  geom_text(aes(label = paste(round(bio1, 2), 'C'), x = 1, y = 15), color = 'black', data = neonsitedata %>% filter(siteID %in% good_sites) %>% mutate(siteID = factor(siteID, levels=sites_temporder$siteID))) +
  theme_john + ggtitle('Mammal Weights ordered by temperature')



ggplot(filter(mam_capture_sitemerge, year == 2015, Pineda_Main_food == 'Herbivore'), aes(x = log10(weight), group = taxonID)) +
  geom_density(adjust = 2) + facet_wrap(~ siteID) +
  theme_minimal() + ggtitle('Herbivores')

ggplot(filter(mam_capture_sitemerge, year == 2015, Pineda_Main_food == 'Generalist'), aes(x = log10(weight), group = taxonID)) +
  geom_density(adjust = 1) + facet_wrap(~ siteID) +
  theme_minimal()+ ggtitle('Generalists')

ggplot(filter(mam_capture_sitemerge, year == 2015, Pineda_Main_food == 'Insectivore'), aes(x = log10(weight), group = taxonID)) +
  geom_density(adjust = 2) + facet_wrap(~ siteID) +
  theme_minimal() + ggtitle('Insectivores')

# Calculate means, variances, and distances -------------------------------

# Mean and variance of the log-transformed weights by species and site.
# Get rid of all with less than 2 individuals

mam_weights <- mam_capture_sitemerge %>%
  filter(year == 2015) %>%
  mutate(logweight = log10(weight)) %>%
  group_by(siteID, taxonID) %>%
  summarize(meanweight = mean(logweight, na.rm=TRUE), sdweight = sd(logweight, na.rm=TRUE), nweight = sum(!is.na(logweight)), diet = Pineda_Main_food[1]) %>%
  filter(nweight > 1)

# Mean of difference to next lowest mammal, and difference to next highest mammal

wtdiff_gran <- mam_weights %>%
  filter(diet %in% c('Generalist','Granivore')) %>%
  arrange(siteID, meanweight) %>%
  group_by(siteID) %>%
  summarize(avgdiff = mean(diff(meanweight)), avgsd = mean(sdweight))
  
ggplot(wtdiff_gran, aes(x = avgdiff, y = avgsd)) + geom_point()  

wtdiff_all <- mam_weights %>%
 # filter(diet %in% c('Generalist','Granivore')) %>%
  arrange(siteID, meanweight) %>%
  group_by(siteID) %>%
  summarize(avgdiff = mean(diff(meanweight)), avgsd = mean(sdweight))

ggplot(wtdiff_all, aes(x = avgdiff, y = avgsd)) + geom_point()  

wtdiff_herb <- mam_weights %>%
  filter(diet %in% c('Generalist','Herbivore')) %>%
  arrange(siteID, meanweight) %>%
  group_by(siteID) %>%
  summarize(avgdiff = mean(diff(meanweight)), avgsd = mean(sdweight))

ggplot(wtdiff_herb, aes(x = avgdiff, y = avgsd)) + geom_point()  