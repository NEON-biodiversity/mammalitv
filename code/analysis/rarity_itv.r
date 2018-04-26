# Exploratory commonness/rarity analysis
# Test predictions made by Umana et al. 2015
# Author: QDR
# Project: NEON ITV
# Created:18 Jan 2017

# Modified 20 Apr 2018: Specify paths to data on HPCC

data_path <- '/mnt/research/neon'

# Load and clean data -----------------------------------------------------

# Use the same pre-digested dataset that is used for o-statistic analysis
# However, use 2016 since there is no reason to exclude it at this point.

library(dplyr)
library(ggplot2)
library(lubridate)
library(reshape2)

load(file.path(data_path, 'final_data/allorganismal_latest.r'))
load(file.path(data_path, 'final_data/allsiteplot_latest.r'))

mammalTraits <- read.csv(file.path(data_path, 'external_data/final_external_data/NEON_miscmammaltraits.csv'))

mam_use <- mutate(mam_capture, 
                      individualandtag = pmin(as.character(individualID), as.character(tagID), na.rm=TRUE),
                      year = year(date)) %>%
  filter(order == 'Rodentia', taxonProtocolCategory == 'target')

# Include guilds in case needed
mammalGuilds <- mammalTraits %>% 
  mutate(scientificName = paste(Genus, Species)) %>%
  dplyr::select(scientificName, Pineda_Main_food)

mammalGuilds <- rbind(mammalGuilds, data.frame(scientificName=c("Dipodomys sp.", "Glaucomys sp.", "Microtus sp.", "Neotoma sp.", 
                                                                "Perognathus sp.", "Peromyscus sp.", "Reithrodontomys sp.", "Sigmodon sp.",
                                                                "Chaetodipus sp.", "Onychomys sp."),
                                               Pineda_Main_food=c('Generalist','Generalist','Herbivore','Herbivore','Generalist','Granivore',
                                                                  'Generalist','Generalist','Granivore','Insectivore')))


mam_use <- left_join(mam_use, mammalGuilds, by = 'scientificName')

# Convert all multiple measured individuals into single measurements
mam_noID <- mam_use[mam_use$individualandtag=='', c('year', 'siteID','plotID', 'taxonID','scientificName','family', 'individualandtag','hindfootLength','earLength','tailLength','totalLength','weight','sex','lifeStage','Pineda_Main_food')]

mam_grp <- filter(mam_use, individualandtag != '') %>%
  group_by(year, siteID, plotID, taxonID, scientificName, family, individualandtag)

mam_byindiv <- summarize_at(mam_grp, vars(hindfootLength, earLength, tailLength, totalLength, weight), median, na.rm=TRUE)
mam_byindiv_class <- do(mam_grp, data.frame(sex = .$sex[1], lifeStage = names(sort(table(.$lifeStage),decreasing=TRUE))[1], Pineda_Main_food=.$Pineda_Main_food[1]))

mam_byindiv <- cbind(as.data.frame(mam_byindiv), as.data.frame(mam_byindiv_class[,c('sex', 'lifeStage','Pineda_Main_food')]))                  
mam_byindiv <-rbind(mam_byindiv, mam_noID)

# If needed, covariates can be joined later.
#mam_merge <- left_join(mam_byindiv, filter(neonplotdata, subtype == 'mammalGrid'))



# Calculate predictor and response ----------------------------------------



# Calculate relative abundance and ITV for all the species.
# Pool by site (not plot)

# for now, include both sexes and life stages.
# Remove all individuals not identified to species.
# Use all years' data

do_rar_itv <- function(dat) {
  dat %>% 
    group_by(taxonID, scientificName, family, Pineda_Main_food) %>%
    do(data.frame(abund = length(.$logweight)/nrow(dat), itv = sd(.$logweight, na.rm=TRUE)/mean(.$logweight, na.rm=TRUE)))
}

# Exclude sites with low trapping success
# Sites with less than 200 records include the Alaska sites and the swamp sites.
nrecord_bysite <- table(mam_byindiv$siteID)
excl_sites <- names(which(nrecord_bysite < 200))

mam_rarity <- mam_byindiv %>%
  filter(!grepl('sp\\.', scientificName), !siteID %in% excl_sites) %>%
  mutate(logweight = log10(weight)) %>%
  group_by(siteID) %>%
  do(do_rar_itv(.)) %>%
  rename(guild = Pineda_Main_food)
    

# Plot results ------------------------------------------------------------



# Plot all together
ggplot(mam_rarity) +
  geom_point(aes(x = abund, y = itv, color = siteID, shape = guild)) +
  theme_minimal() +
  labs(x = 'Relative abundance', y = 'CV of log10 body mass')

# Plot them by guild
ggplot(mam_rarity) +
  geom_point(aes(x = abund, y = itv, color = siteID)) +
  facet_wrap(~ guild) +
  theme_minimal() +
  labs(x = 'Relative abundance', y = 'CV of log10 body mass')

# Plot them by site
ggplot(mam_rarity) +
  geom_point(aes(x = abund, y = itv, shape = guild)) +
  facet_wrap(~ siteID) +
  theme_minimal() +
  labs(x = 'Relative abundance', y = 'CV of log10 body mass')


# Hypothesis testing ------------------------------------------------------



# Run a simple model on them
# Site is a random effect.
# Must transform the relative abundances and the ITVs. Logit transformation on both since they are both proportions?

library(lme4)
rarity_lmm <- lmer(itv_trans ~ abund_trans + (1|siteID), data = mam_rarity %>% mutate(itv_trans = log(itv), abund_trans = qlogis(abund)))

summary(rarity_lmm)
confint(rarity_lmm, method = 'boot', nsim = 999, parm = 'abund_trans') # Bootstrap confidence interval of fixed effect coefficient. Not significant at 95%
confint(rarity_lmm, method = 'boot', nsim = 999, parm = 'abund_trans', level = 0.90) # Also not quite significant at 90%, though it is trending up.

library(MuMIn)
r.squaredGLMM(rarity_lmm) # Only 1% of the variation is explained by the fixed effect.


# Test prediction about heterogeneity -------------------------------------

# The prediction is that the relationship will flip such that in homogeneous environments, common species will have less ITV (they're common because they've adapted to the optimum).
# In heterogeneous environments, common species will have more ITV (they're common because they can deal with the heterogeneity)

# So the prediction is that as heterogeneity increases, the slope of the itv~abundance relationship will go from negative to positive

# Load dataframe of environmental heterogeneity
# We do not have site data for CLBJ yet, but that might be OK. Would need to go back and redownload neon site data.
heterodf <- read.csv(file.path(data_path, 'MS1_RodentOverlap/R_data/heterogeneity.csv', stringsAsFactors = FALSE))

# Test whether the slope of the itv~abundance relationship is a function of heterogeneity.

itv_abund_slopes <- mam_rarity %>% 
  mutate(itv_trans = log(itv), abund_trans = qlogis(abund)) %>%
  group_by(siteID) %>% 
  do(data.frame(slope = lm(.$itv_trans ~ .$abund_trans)$coef[2])) %>%
  left_join(heterodf)

ggplot(itv_abund_slopes, aes(x = pc1_productivityheterogeneity, y = slope)) + 
  geom_point() + theme_minimal() + labs(x = 'environmental heterogeneity', y = 'slope of relationship')

ggplot(itv_abund_slopes, aes(x = pc2_topographyheterogeneity, y = slope)) + 
  geom_point() + theme_minimal() + labs(x = 'topographic heterogeneity', y = 'slope of relationship')

slope_lm <- lm(slope ~ pc1_productivityheterogeneity + pc2_topographyheterogeneity, data = itv_abund_slopes)

summary(slope_lm) # The prediction is not supported for mammal body size.
