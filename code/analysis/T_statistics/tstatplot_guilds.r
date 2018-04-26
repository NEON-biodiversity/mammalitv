# Plot mammal t-stats, separated out by guild
# Author: QDR
# Project: NEON ITV
# Created: 22 Sep 2016
# Last modified: 26 Sep 2016

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
load('C:/Users/Q/Dropbox/neon/code/tstats_2015_gg_iucn.r')
load('C:/Users/Q/Dropbox/neon/code/tstats_2015_gh_iucn.r')
load('C:/Users/Q/Dropbox/neon/code/tstats_2015_ggh_iucn.r')
load('C:/Users/Q/Dropbox/neon/code/mammalPDbyplotobject.r')

mammalTax <- read.csv('C:/Users/Q/Dropbox/neon/protocols/taxonomy/NEON_mam_taxonomy.csv')
mammalTraits <- read.csv('C:/Users/Q/Dropbox/neon/data/external_datasets/NEON_miscmammaltraits.csv')


mam_capture <- mutate(mam_capture, 
                      individualandtag = pmin(as.character(individualID), as.character(tagID), na.rm=TRUE),
                      year = year(date)) %>%
  left_join(mammalTax[,c('taxonID','order')]) %>% # Added 21 Sep. Gets rid of everything but rodents.
  filter(order == 'Rodentia')

mammalGuilds <- mammalTraits %>% 
  mutate(scientificName = paste(Genus, Species)) %>%
  select(scientificName, Pineda_Main_food)

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

tmgg <- tmgg %>% left_join(chao1gg15) %>% left_join(mammalPD)
tmgh <- tmgh %>% left_join(chao1gg15) %>% left_join(mammalPD)
tmggh <- tmggh %>% left_join(chao1gg15) %>% left_join(mammalPD)



# Create scatterplots -----------------------------------------------------

# Remove outlier 
tmgg <- filter(tmgg, !(stat == 'T_IP.IC' & t < -10))
tmgh <- filter(tmgh, !(stat == 'T_IP.IC' & t < -10))
tmggh <- filter(tmggh, !(stat == 'T_IP.IC' & t < -10))

pttempgg <- ggplot(subset(tmgg, trait == 'logweight'), aes(x=bio1)) + 
  facet_wrap( ~ stat, scales = 'free') +
  geom_segment(aes(y = ci_min, yend = ci_max, xend=bio1), alpha = 0.5, size = 1.5, color = 'skyblue') +
  geom_point(aes(y = t), size = 1.5) +
  labs(y = 'T-statistic', x = parse(text=bioclimnames[1])) +
  theme_minimal() + theme(panel.border = element_rect(fill='transparent')) +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'midnightblue', size=0.9) +
  qSubtitle('T-statistics for 2015 NEON mammals versus MAT', 'Generalists + Granivores')

ptmpdgg <- ggplot(subset(tmgg, trait == 'logweight'), aes(x=mpd.obs.z2015)) + 
  facet_grid(~ stat, scales = 'free') +
  geom_segment(aes(y = ci_min, yend = ci_max, xend=mpd.obs.z2015), alpha = 0.5, size = 1.5, color = 'skyblue') +
  geom_point(aes(y = t), size = 1.5) +
  labs(y = 'T-statistic', x = 'SES of mean pairwise distance') +
  theme_minimal() + theme(panel.border = element_rect(fill='transparent')) +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'midnightblue', size=0.9) +
  qSubtitle('T-statistics for 2015 NEON mammals versus MPD', 'Generalists + Granivores')

ptchaogg <- ggplot(subset(tmgg, trait == 'logweight'), aes(x=chao1)) + 
  facet_grid(~ stat, scales = 'free') +
  geom_segment(aes(y = ci_min, yend = ci_max, xend=chao1), alpha = 0.5, size = 1.5, color = 'skyblue') +
  geom_point(aes(y = t), size = 1.5) +
  labs(y = 'T-statistic', x = 'SES of mean pairwise distance') +
  theme_minimal() + theme(panel.border = element_rect(fill='transparent')) +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'midnightblue', size=0.9) +
  qSubtitle('T-statistics for 2015 NEON mammals versus richness', 'Generalists + Granivores')

pttempgh <- ggplot(subset(tmgh, trait == 'logweight'), aes(x=bio1)) + 
  facet_wrap( ~ stat, scales = 'free') +
  geom_segment(aes(y = ci_min, yend = ci_max, xend=bio1), alpha = 0.5, size = 1.5, color = 'skyblue') +
  geom_point(aes(y = t), size = 1.5) +
  labs(y = 'T-statistic', x = parse(text=bioclimnames[1])) +
  theme_minimal() + theme(panel.border = element_rect(fill='transparent')) +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'midnightblue', size=0.9) +
  qSubtitle('T-statistics for 2015 NEON mammals versus MAT', 'Generalists + Herbivores')

ptmpdgh <- ggplot(subset(tmgh, trait == 'logweight'), aes(x=mpd.obs.z2015)) + 
  facet_grid(~ stat, scales = 'free') +
  geom_segment(aes(y = ci_min, yend = ci_max, xend=mpd.obs.z2015), alpha = 0.5, size = 1.5, color = 'skyblue') +
  geom_point(aes(y = t), size = 1.5) +
  labs(y = 'T-statistic', x = 'SES of mean pairwise distance') +
  theme_minimal() + theme(panel.border = element_rect(fill='transparent')) +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'midnightblue', size=0.9) +
  qSubtitle('T-statistics for 2015 NEON mammals versus MPD', 'Generalists + Herbivores')

ptchaogh <- ggplot(subset(tmgh, trait == 'logweight'), aes(x=chao1)) + 
  facet_grid(~ stat, scales = 'free') +
  geom_segment(aes(y = ci_min, yend = ci_max, xend=chao1), alpha = 0.5, size = 1.5, color = 'skyblue') +
  geom_point(aes(y = t), size = 1.5) +
  labs(y = 'T-statistic', x = 'SES of mean pairwise distance') +
  theme_minimal() + theme(panel.border = element_rect(fill='transparent')) +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'midnightblue', size=0.9) +
  qSubtitle('T-statistics for 2015 NEON mammals versus richness', 'Generalists + Herbivores')

pttempggh <- ggplot(subset(tmggh, trait == 'logweight'), aes(x=bio1)) + 
  facet_wrap( ~ stat, scales = 'free') +
  geom_segment(aes(y = ci_min, yend = ci_max, xend=bio1), alpha = 0.5, size = 1.5, color = 'skyblue') +
  geom_point(aes(y = t), size = 1.5) +
  labs(y = 'T-statistic', x = parse(text=bioclimnames[1])) +
  theme_minimal() + theme(panel.border = element_rect(fill='transparent')) +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'midnightblue', size=0.9) +
  qSubtitle('T-statistics for 2015 NEON mammals versus MAT', 'Generalists + Granivores + Herbivores')

ptmpdggh <- ggplot(subset(tmggh, trait == 'logweight'), aes(x=mpd.obs.z2015)) + 
  facet_grid(~ stat, scales = 'free') +
  geom_segment(aes(y = ci_min, yend = ci_max, xend=mpd.obs.z2015), alpha = 0.5, size = 1.5, color = 'skyblue') +
  geom_point(aes(y = t), size = 1.5) +
  labs(y = 'T-statistic', x = 'SES of mean pairwise distance') +
  theme_minimal() + theme(panel.border = element_rect(fill='transparent')) +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'midnightblue', size=0.9) +
  qSubtitle('T-statistics for 2015 NEON mammals versus MPD', 'Generalists + Granivores + Herbivores')

ptchaoggh <- ggplot(subset(tmggh, trait == 'logweight'), aes(x=chao1)) + 
  facet_grid(~ stat, scales = 'free') +
  geom_segment(aes(y = ci_min, yend = ci_max, xend=chao1), alpha = 0.5, size = 1.5, color = 'skyblue') +
  geom_point(aes(y = t), size = 1.5) +
  labs(y = 'T-statistic', x = 'SES of mean pairwise distance') +
  theme_minimal() + theme(panel.border = element_rect(fill='transparent')) +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'midnightblue', size=0.9) +
  qSubtitle('T-statistics for 2015 NEON mammals versus richness', 'Generalists + Granivores + Herbivores')

pdf('figs/pdf/tstat_guild_scatterplots.pdf', height = 6, width = 8)
pttempgg
ptmpdgg
ptchaogg

pttempgh
ptmpdgh
ptchaogh

pttempggh
ptmpdggh
ptchaoggh
dev.off()


# Scatterplots by site mean, min, and max.

tmggsumm <- tmgg %>%
  group_by(trait, stat, siteID) %>%
  summarize_at(c('t','bio1','mpd.obs.z2015','chao1','ruggedness'), funs(mean, min, max))

pttempggsite <- ggplot(subset(tmggsumm, trait == 'logweight'), aes(x=bio1_mean)) + 
  facet_wrap( ~ stat, scales = 'free') +
 # geom_segment(aes(y = ci_min, yend = ci_max, xend=bio1), alpha = 0.5, size = 1.5, color = 'skyblue') +
  geom_errorbar(aes(ymin = t_min, ymax = t_max)) +
  geom_point(aes(y = t_mean)) +
  labs(y = 'T-statistic', x = parse(text=bioclimnames[1])) +
  theme_minimal() + theme(panel.border = element_rect(fill='transparent')) +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'blue', size=0.9) +
  qSubtitle('T-statistics for 2015 NEON mammals versus MAT', 'Generalists + Granivores, by site')

ptmpdggsite <- ggplot(subset(tmggsumm, trait == 'logweight'), aes(x=mpd.obs.z2015_mean)) + 
  facet_wrap( ~ stat, scales = 'free') +
  # geom_segment(aes(y = ci_min, yend = ci_max, xend=bio1), alpha = 0.5, size = 1.5, color = 'skyblue') +
  geom_errorbar(aes(ymin = t_min, ymax = t_max)) +
  geom_point(aes(y = t_mean)) +
  labs(y = 'T-statistic', x = 'SES of mean pairwise distance') +
  theme_minimal() + theme(panel.border = element_rect(fill='transparent')) +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'blue', size=0.9) +
  qSubtitle('T-statistics for 2015 NEON mammals versus MPD', 'Generalists + Granivores, by site')

# Run regression on log weights -------------------------------------------

dat_tipic_gg <- tmgg %>% filter(trait == 'logweight', stat == 'T_IP.IC') %>%
  select(trait, t, siteID, plotID, bio1, bio4, cv_bio1, bio12, bio15, cv_bio12, chao1, mpd.obs.z2015, mntd.obs.z2015, ruggedness) %>%
  filter(complete.cases(.)) %>%
  rename(mpd_z = mpd.obs.z2015, mntd_z = mntd.obs.z2015) %>%
  mutate_each(funs((. - mean(.))/sd(.)), bio1:ruggedness) 

dat_all_gg <- tmgg %>% filter(trait == 'logweight') %>%
  select(trait, stat, t, siteID, plotID, bio1, bio4, cv_bio1, bio12, bio15, cv_bio12, chao1, mpd.obs.z2015, mntd.obs.z2015, ruggedness) %>%
  filter(complete.cases(.)) %>%
  rename(mpd_z = mpd.obs.z2015, mntd_z = mntd.obs.z2015) %>%
  mutate_each(funs((. - mean(.))/sd(.)), bio1:ruggedness) 

dat_all_gh <- tmgh %>% filter(trait == 'logweight') %>%
  select(trait, stat, t, siteID, plotID, bio1, bio4, cv_bio1, bio12, bio15, cv_bio12, chao1, mpd.obs.z2015, mntd.obs.z2015, ruggedness) %>%
  filter(complete.cases(.)) %>%
  rename(mpd_z = mpd.obs.z2015, mntd_z = mntd.obs.z2015) %>%
  mutate_each(funs((. - mean(.))/sd(.)), bio1:ruggedness) 

dat_all_ggh <- tmggh %>% filter(trait == 'logweight') %>%
  select(trait, stat, t, siteID, plotID, bio1, bio4, cv_bio1, bio12, bio15, cv_bio12, chao1, mpd.obs.z2015, mntd.obs.z2015, ruggedness) %>%
  filter(complete.cases(.)) %>%
  rename(mpd_z = mpd.obs.z2015, mntd_z = mntd.obs.z2015) %>%
  mutate_each(funs((. - mean(.))/sd(.)), bio1:ruggedness) 

# .... add the other stuff here.

library(brms)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Bayesian

brm_tipic_gg <- brm(t ~ bio1 + chao1 + mpd_z + ruggedness + (1|siteID), data = subset(dat_all_gg, stat == 'T_IP.IC'), chains=4, cores=4, iter = 10000, warmup = 5000)
brm_ticir_gg <- brm(t ~ bio1 + chao1 + mpd_z + ruggedness + (1|siteID), data = subset(dat_all_gg, stat == 'T_IC.IR'), chains=4, cores=4, iter = 10000, warmup = 5000)

sum_tipic_gg <- summary(brm_tipic_gg)

# Frequentist

library(lme4)
library(MuMIn)

glm_tipic_gg <- lmer(t ~ bio1 + chao1 + mpd_z + ruggedness + (1|siteID), 
                     data = dat_all_gg, subset = stat == 'T_IP.IC', na.action = 'na.pass', REML=FALSE) # Very similar to bayesian
summary(glm_tipic_gg)
confint(glm_tipic_gg, method='boot', nsim = 999)
dredge_tipic_gg <- dredge(glm_tipic_gg)
r.squaredGLMM(glm_tipic_gg)

# Extract best models below preset threshold
tipic_ggbest <- get.models(dredge_tipic_gg, subset = delta < 1)

# Run model averaging on all models below threshold
tipic_ggavg <- model.avg(tipic_ggbest)

# Show which parameters are in the best models
plot(subset(dredge_tipic_gg, delta < 1))

# Show summary of parameter estimates in single best model only
summary(tipic_ggbest[[1]]) 

# Get bootstrap CI of parameters
confint(tipic_ggbest[[1]], method = 'boot', nsim = 999)

# Find R-squared of the fixed and random effects
r.squaredGLMM(tipic_ggbest[[1]])


### Run all models together.

data_list <- list(dat_all_gg, dat_all_gh, dat_all_ggh)
stat_list <- list('T_IP.IC', 'T_IC.IR', 'T_PC.PR')
fit_list <- list()
confint_list <- list()
coeff_list <- list()
rsq_list <- list()

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
r2data <- data.frame(r2data, community = rep(c('Gen + Gran', 'Gen + Herb', 'Gen + Gran + Herb'), each = 3), stat = c('T_IP.IC','T_IC.IR','T_PC.PR'))

library(reshape2)
cimins <- melt(cimins, varnames = c('predictor','model'), value.name = 'ci_min')
cimaxes <- melt(cimaxes, varnames = c('predictor','model'), value.name = 'ci_max')
coeffs <- melt(coeffs, varnames = c('predictor','model'), value.name = 'coefficient')

fits_dat <- coeffs %>% left_join(cimins) %>% left_join(cimaxes) %>%
  mutate(community = rep(c('Gen + Gran', 'Gen + Herb', 'Gen + Gran + Herb'), each = 12),
         stat = rep(rep(c('T_IP.IC','T_IC.IR','T_PC.PR'), each = 4), times=3))

library(ggplot2)

ggplot(left_join(fits_dat, r2data), aes(x = predictor, y = coefficient, ymin = ci_min, ymax = ci_max)) +
  facet_grid(community ~ stat) + 
  geom_hline(yintercept = 0, color = 'blue', linetype = 'dotted') +
  geom_errorbar(width = 0.1) +
  geom_point() +
  geom_text(aes(x = 4.5, y = -1, label = round(R2m,2))) +
  coord_flip() +
  theme_bw()

ggsave('figs/png/tstat_guild_coefficients.png', height=9, width=9, dpi=400)