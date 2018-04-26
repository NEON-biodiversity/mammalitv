# Fit Bayesian regression model to T-stats (NEON mammals) using packages brms or rstanarm
# Author: QDR
# Project: NEON ITV
# Created: 01 Sep 2016
# Last modified: 06 Sep 2016

# Modified 06 Sep 2016: Added improved spatial stats.

# Check out https://thinkinator.com/2016/01/12/r-users-will-now-inevitably-become-bayesians/

# The brms and rstanarm packages are front-ends to stan that allow similar syntax to glm in R.

# Load analysis packages
#library(rstanarm)

# Load data

library(dplyr)
library(lubridate)
library(cati)
source('~/GitHub/NEON/code/analysis/tstat2longform.r')

load('~/GitHub/NEON/allorganismal_latest.r')
load('~/GitHub/NEON/allsiteplot_latest.r')
load('C:/Users/Q/Dropbox/neon/code/tstats_mam_plots15_iucnpool.r')
load('C:/Users/Q/Dropbox/neon/code/mammalPDbyplotobject.r')

spatialstat <- read.csv('C:/Users/Q/Dropbox/neon/data/external_datasets/NEON_spatial_stats_highres.csv')
neonplotdata <- cbind(neonplotdata, spatialstat) %>% rename(ruggedness = tri)

tmiucn <- tstat2longform_ses(tstats_mam_plots15_iucnpool, bysite = FALSE) %>% 
  left_join(filter(neonplotdata, subtype == 'mammalGrid'))

# Mammals have abundance so we use the Chao1 estimator.

chao <- function(x) {
  xcomm <- table(x$taxonID)
  S_obs <- length(xcomm)
  f1 <- sum(xcomm == 1)
  f2 <- sum(xcomm == 2)
  return(data.frame(chao1 = S_obs + (f1 * (f1 - 1)) / (2 * (f2 + 1))))
}

chao1mammal15 <- mam_capture %>% filter(year(date) == 2015) %>%
  group_by(plotID) %>% do(chao(.))

tmiucn <- left_join(tmiucn, chao1mammal15)
tmiucn <- left_join(tmiucn, mammalPD)


# Create data objects.

dat_tipic_wt <- tmiucn %>% filter(trait == 'weight', stat == 'T_IP.IC') %>%
  select(t, siteID, bio1, bio4, cv_bio1, bio12, bio15, cv_bio12, chao1, mpd.obs.z2015, mntd.obs.z2015, ruggedness) %>%
  filter(complete.cases(.)) %>%
  rename(mpd_z = mpd.obs.z2015, mntd_z = mntd.obs.z2015) %>%
  #mutate(siteID = as.numeric(factor(siteID))) %>%
  mutate_each(funs((. - mean(.))/sd(.)), bio1:ruggedness) 

dat_tipic_hf <- tmiucn %>% filter(trait == 'hindfootLength', stat == 'T_IP.IC') %>%
  select(t, siteID, bio1, bio4, cv_bio1, bio12, bio15, cv_bio12, chao1, mpd.obs.z2015, mntd.obs.z2015, ruggedness) %>%
  filter(complete.cases(.)) %>%
  rename(mpd_z = mpd.obs.z2015, mntd_z = mntd.obs.z2015) %>%
  #mutate(siteID = as.numeric(factor(siteID))) %>%
  mutate_each(funs((. - mean(.))/sd(.)), bio1:ruggedness) 

dat_ticir_wt <- tmiucn %>% filter(trait == 'weight', stat == 'T_IC.IR') %>%
  select(t, siteID, bio1, bio4, cv_bio1, bio12, bio15, cv_bio12, chao1, mpd.obs.z2015, mntd.obs.z2015, ruggedness) %>%
  filter(complete.cases(.)) %>%
  rename(mpd_z = mpd.obs.z2015, mntd_z = mntd.obs.z2015) %>%
  #mutate(siteID = as.numeric(factor(siteID))) %>%
  mutate_each(funs((. - mean(.))/sd(.)), bio1:ruggedness) 

dat_ticir_hf <- tmiucn %>% filter(trait == 'hindfootLength', stat == 'T_IC.IR') %>%
  select(t, siteID, bio1, bio4, cv_bio1, bio12, bio15, cv_bio12, chao1, mpd.obs.z2015, mntd.obs.z2015, ruggedness) %>%
  filter(complete.cases(.)) %>%
  rename(mpd_z = mpd.obs.z2015, mntd_z = mntd.obs.z2015) %>%
  #mutate(siteID = as.numeric(factor(siteID))) %>%
  mutate_each(funs((. - mean(.))/sd(.)), bio1:ruggedness) 

dat_tpcpr_wt <- tmiucn %>% filter(trait == 'weight', stat == 'T_PC.PR') %>%
  select(t, siteID, bio1, bio4, cv_bio1, bio12, bio15, cv_bio12, chao1, mpd.obs.z2015, mntd.obs.z2015, ruggedness) %>%
  filter(complete.cases(.)) %>%
  rename(mpd_z = mpd.obs.z2015, mntd_z = mntd.obs.z2015) %>%
  #mutate(siteID = as.numeric(factor(siteID))) %>%
  mutate_each(funs((. - mean(.))/sd(.)), bio1:ruggedness) 

dat_tpcpr_hf <- tmiucn %>% filter(trait == 'hindfootLength', stat == 'T_PC.PR') %>%
  select(t, siteID, bio1, bio4, cv_bio1, bio12, bio15, cv_bio12, chao1, mpd.obs.z2015, mntd.obs.z2015, ruggedness) %>%
  filter(complete.cases(.)) %>%
  rename(mpd_z = mpd.obs.z2015, mntd_z = mntd.obs.z2015) %>%
  #mutate(siteID = as.numeric(factor(siteID))) %>%
  mutate_each(funs((. - mean(.))/sd(.)), bio1:ruggedness) 

# Multiple response variable regression
dat_tipic_all <- tmiucn %>% filter(trait %in% c('weight', 'hindfootLength'), stat == 'T_IP.IC') %>%
  select(trait, t, siteID, plotID, bio1, bio4, cv_bio1, bio12, bio15, cv_bio12, chao1, mpd.obs.z2015, mntd.obs.z2015, ruggedness) %>%
  filter(complete.cases(.)) %>%
  rename(mpd_z = mpd.obs.z2015, mntd_z = mntd.obs.z2015) %>%
  #mutate(siteID = as.numeric(factor(siteID))) %>%
  mutate_each(funs((. - mean(.))/sd(.)), bio1:ruggedness) %>%
  split(.$trait, drop =TRUE)
dat_tipic_all <- left_join(rename(dat_tipic_all[[1]], hindfootLength_t = t)[,-1], rename(dat_tipic_all[[2]], weight_t = t)[,2:4])  

dat_ticir_all <- tmiucn %>% filter(trait %in% c('weight', 'hindfootLength'), stat == 'T_IC.IR') %>%
  select(trait, t, siteID, plotID, bio1, bio4, cv_bio1, bio12, bio15, cv_bio12, chao1, mpd.obs.z2015, mntd.obs.z2015, ruggedness) %>%
  filter(complete.cases(.)) %>%
  rename(mpd_z = mpd.obs.z2015, mntd_z = mntd.obs.z2015) %>%
  #mutate(siteID = as.numeric(factor(siteID))) %>%
  mutate_each(funs((. - mean(.))/sd(.)), bio1:ruggedness) %>%
  split(.$trait, drop =TRUE)
dat_ticir_all <- left_join(rename(dat_ticir_all[[1]], hindfootLength_t = t)[,-1], rename(dat_ticir_all[[2]], weight_t = t)[,2:4])

dat_tpcpr_all <- tmiucn %>% filter(trait %in% c('weight', 'hindfootLength'), stat == 'T_PC.PR') %>%
  select(trait, t, siteID, plotID, bio1, bio4, cv_bio1, bio12, bio15, cv_bio12, chao1, mpd.obs.z2015, mntd.obs.z2015, ruggedness) %>%
  filter(complete.cases(.)) %>%
  rename(mpd_z = mpd.obs.z2015, mntd_z = mntd.obs.z2015) %>%
  #mutate(siteID = as.numeric(factor(siteID))) %>%
  mutate_each(funs((. - mean(.))/sd(.)), bio1:ruggedness) %>%
  split(.$trait, drop =TRUE)
dat_tpcpr_all <- left_join(rename(dat_tpcpr_all[[1]], hindfootLength_t = t)[,-1], rename(dat_tpcpr_all[[2]], weight_t = t)[,2:4])

# "cheating" multivariate regression

# library(lme4)
# tipiclmer <- lmer(t ~ trait:bio1 + trait:mpd_z + (1|siteID), data = dat_tipic_all, na.action = 'na.pass', REML = FALSE)


############## below, attempt to split and run multivariate

# tipicsplit <- split(dat_tipic_all, dat_tipic_all$trait, drop = TRUE)
# dat_tipic_all <- left_join(rename(tipicsplit[[1]], hindfootLength_t = t)[,-1], rename(tipicsplit[[2]], weight_t = t)[,2:4])



# so far, none of the stuff below works.

# Multiple responses
# stan_glmer(cbind(hindfootLength_t, weight_t) ~ bio1 + bio4 + (1|siteID), data = dat_tipic_all, prior = normal())
# 
# # Multivariate GLMM
# # the us(trait):siteID is the most conservative and highly parameterized random effects specification (can be simplified by changing assumptions). It allows for heterogeneous variance for each trait, and covariance among traits. "trait" is a reserved word in the mcmcglmm package referring to a response variable.
# # The us(trait):units is the fully parameterized residual covariance matrix (can be simplified)
# # Default priors are used here
# 
# prior = list(R = list(V = diag(2)/3, n = 2),
#              G = list(lapply(1:24, function(x) list(V=diag(2)/3, n=2))))
# 
# prior = list(R = list(V = diag(2)/3, n = 2),
#              G = list(G1 = list(V = diag(2)/3, n = 2),
#                       G2 = list(V = diag(2)/3, n = 2)))
# 
# tipicglmm <- MCMCglmm(cbind(hindfootLength_t, weight_t) ~ trait:bio1 + trait:chao1 - 1, 
#                       random = ~ us(trait):siteID,
#                       rcov = ~ us(trait):units,
#                       prior = prior,
#                       family = c('gaussian', 'gaussian'),
#                       data = select(dat_tipic_all, hindfootLength_t, weight_t, bio1, chao1, siteID), 
#                       nitt = 10000,
#                       burnin = 1000,
#                       thin = 1)

library(brms)
brm(weight_t ~ bio1 + chao1 + (1|siteID), data = dat_tipic_all)
brm(cbind(weight_t, hindfootLength_t) ~ bio1 + chao1 + (1|siteID), data = dat_tipic_all) # Not correct, but it runs.

rstan_options(auto_write=TRUE)
option(mc.cores=2)

brm1 <- brm(cbind(weight_t, hindfootLength_t) ~ 0 + trait + trait:(bio1 + bio4 + cv_bio1 + bio12 + bio15 + cv_bio12 + chao1 + mpd_z + mntd_z + ruggedness) + (1|siteID), data = dat_tipic_all, chains=4, cores=4, iter = 10000, warmup = 5000)

# Try to reduce brm model.
waic.full <- WAIC(brm1) # Calculate information criterion for model.

# Get rid of predictors
# Pairs plot: look at collinearity of predictors and get rid of the ones not needed
# eta should be the fixed effect predictors.

# should be faster to write to file
png('C:/Users/Q/Desktop/brmpairs.png', height=12, width=12, units = 'in', res = 300)
pairs(brm1, pars = 'b')
dev.off()

# Reduced model
brmreduced_tipic <- brm(cbind(weight_t, hindfootLength_t) ~ 0 + trait:(bio1 + chao1 + mpd_z + ruggedness) + (1|siteID), data = dat_tipic_all, chains=4, cores=4, iter = 10000, warmup = 5000)


# png('~/figs/brms/brmreduced_tipic.png', height=12, width=12, units='in', res=300)
# pairs(brmreduced_tipic, pars='b')
# dev.off()

brmreduced_ticir <- brm(cbind(weight_t, hindfootLength_t) ~ 0 + trait:(bio1 + chao1 + mpd_z + ruggedness) + (1|siteID), data = dat_ticir_all, chains=4, cores=4, iter = 10000, warmup = 5000)
brmreduced_tpcpr <- brm(cbind(weight_t, hindfootLength_t) ~ 0 + trait:(bio1 + chao1 + mpd_z + ruggedness) + (1|siteID), data = dat_tpcpr_all, chains=4, cores=4, iter = 10000, warmup = 5000)

sum_tipic <- summary(brmreduced_tipic)
sum_ticir <- summary(brmreduced_ticir)
sum_tpcpr <- summary(brmreduced_tpcpr)

save(brmreduced_tipic, brmreduced_ticir, brmreduced_tpcpr, file = 'C:/Users/Q/Dropbox/neon/code/brmreducedobjects.r')

# Make plots
summ2plotdat <- function(summ) {
  res <- summ$fixed
  params <- strsplit(dimnames(res)[[1]], ':')
  res <- data.frame(predictor = sapply(params, '[', 2),
             trait = sapply(params, '[', 1),
             res)
  dplyr::rename(res, ci_lower = l.95..CI, ci_upper = u.95..CI)
}

ggplot(summ2plotdat(sum_tipic), aes(x = predictor, y = Estimate, ymin = ci_lower, ymax = ci_upper)) +
  geom_pointrange() + theme_minimal() + 
  facet_grid(. ~ trait) +
  coord_flip()

pdat <- rbind(transform(summ2plotdat(sum_tipic), stat='Tipic'), transform(summ2plotdat(sum_ticir), stat='Ticir'))
pdat <- transform(pdat, not_zero = (ci_lower<0 & ci_upper<0) | (ci_lower>0 & ci_upper>0))

lab1 <- labeller(stat = c(Tipic = 'Individuals within communities', Ticir = 'Communities within regions'), trait = c(traithindfootLength_t = 'Hind foot length', traitweight_t = 'Body weight'))

library(extrafont)

pcoef <- ggplot(pdat, aes(x = predictor, y = Estimate, ymin = ci_lower, ymax = ci_upper)) +
  geom_hline(yintercept = 0, color = 'blue', linetype = 'dotted') +
  geom_errorbar(width = 0.1) +
  geom_point(aes(color = not_zero, size = not_zero)) +
  facet_grid(stat ~ trait, labeller = lab1) +
  coord_flip() +
  scale_color_manual(values = c('black', 'indianred')) +
  scale_size_manual(values = c(2,3)) +
  scale_x_discrete(labels = c('Temperature', 'Species\nrichness', 'Phylogenetic\ndiversity', 'Topographic\nheterogeneity')) +
  theme_minimal() + theme(legend.position = 'none', panel.border = element_rect(fill = 'transparent'), text = element_text(family = 'Century Gothic'), axis.title.y = element_blank())

ggsave('C:/Users/Q/Dropbox/presentations/hanoverSep2016/neonfigs/multivariate_coeff_plot.png', height = 6, width = 7, dpi = 400)  
