# Fit Bayesian regression model to T-stats (NEON mammals)
# Author: QDR
# Project: NEON ITV
# Created: 31 Aug 2016
# Last modified:

# Load stan
library(rstan)
rstan_options(auto_write = TRUE) # This line not needed if run on cluster
options(mc.cores = 4) # Modify this as needed

# Compile model
tstat_model <- stan_model('code/analysis/tstat_multiple_reg.stan')

# Define data objects (for each trait/stat combination)

# Load data

library(dplyr)
library(lubridate)
library(cati)
source('~/GitHub/NEON/code/analysis/tstat2longform.r')

load('~/GitHub/NEON/allorganismal_latest.r')
load('~/GitHub/NEON/allsiteplot_latest.r')
load('C:/Users/Q/Dropbox/neon/code/tstats_mam_plots15_iucnpool.r')
load('C:/Users/Q/Dropbox/neon/code/mammalPDbyplotobject.r')

spatialstat <- read.csv('C:/Users/Q/Dropbox/neon/data/external_datasets/NEON_spatial_stats.csv')
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
  mutate(siteID = as.numeric(factor(siteID))) %>%
  mutate_each(funs((. - mean(.))/sd(.)), bio1:ruggedness) %>%
  c(N = nrow(.), Nsites = max(.$siteID), Nparams = 10)

dat_tipic_hf <- tmiucn %>% filter(trait == 'hindfootLength', stat == 'T_IP.IC') %>%
  select(t, siteID, bio1, bio4, cv_bio1, bio12, bio15, cv_bio12, chao1, mpd.obs.z2015, mntd.obs.z2015, ruggedness) %>%
  filter(complete.cases(.)) %>%
  rename(mpd_z = mpd.obs.z2015, mntd_z = mntd.obs.z2015) %>%
  mutate(siteID = as.numeric(factor(siteID))) %>%
  mutate_each(funs((. - mean(.))/sd(.)), bio1:ruggedness) %>%
  c(N = nrow(.), Nsites = max(.$siteID), Nparams = 10)

dat_ticir_wt <- tmiucn %>% filter(trait == 'weight', stat == 'T_IC.IR') %>%
  select(t, siteID, bio1, bio4, cv_bio1, bio12, bio15, cv_bio12, chao1, mpd.obs.z2015, mntd.obs.z2015, ruggedness) %>%
  filter(complete.cases(.)) %>%
  rename(mpd_z = mpd.obs.z2015, mntd_z = mntd.obs.z2015) %>%
  mutate(siteID = as.numeric(factor(siteID))) %>%
  mutate_each(funs((. - mean(.))/sd(.)), bio1:ruggedness) %>%
  c(N = nrow(.), Nsites = max(.$siteID), Nparams = 10)

dat_ticir_hf <- tmiucn %>% filter(trait == 'hindfootLength', stat == 'T_IC.IR') %>%
  select(t, siteID, bio1, bio4, cv_bio1, bio12, bio15, cv_bio12, chao1, mpd.obs.z2015, mntd.obs.z2015, ruggedness) %>%
  filter(complete.cases(.)) %>%
  rename(mpd_z = mpd.obs.z2015, mntd_z = mntd.obs.z2015) %>%
  mutate(siteID = as.numeric(factor(siteID))) %>%
  mutate_each(funs((. - mean(.))/sd(.)), bio1:ruggedness) %>%
  c(N = nrow(.), Nsites = max(.$siteID), Nparams = 10)

dat_tpcpr_wt <- tmiucn %>% filter(trait == 'weight', stat == 'T_PC.PR') %>%
  select(t, siteID, bio1, bio4, cv_bio1, bio12, bio15, cv_bio12, chao1, mpd.obs.z2015, mntd.obs.z2015, ruggedness) %>%
  filter(complete.cases(.)) %>%
  rename(mpd_z = mpd.obs.z2015, mntd_z = mntd.obs.z2015) %>%
  mutate(siteID = as.numeric(factor(siteID))) %>%
  mutate_each(funs((. - mean(.))/sd(.)), bio1:ruggedness) %>%
  c(N = nrow(.), Nsites = max(.$siteID), Nparams = 10)

dat_tpcpr_hf <- tmiucn %>% filter(trait == 'hindfootLength', stat == 'T_PC.PR') %>%
  select(t, siteID, bio1, bio4, cv_bio1, bio12, bio15, cv_bio12, chao1, mpd.obs.z2015, mntd.obs.z2015, ruggedness) %>%
  filter(complete.cases(.)) %>%
  rename(mpd_z = mpd.obs.z2015, mntd_z = mntd.obs.z2015) %>%
  mutate(siteID = as.numeric(factor(siteID))) %>%
  mutate_each(funs((. - mean(.))/sd(.)), bio1:ruggedness) %>%
  c(N = nrow(.), Nsites = max(.$siteID), Nparams = 10)

# Set sampling options
n_chains <- 4
n_iter <- 9999
n_warmup <- 1000
n_thin <- 1


# Sample model
fit_tipic_wt <- sampling(tstat_model, data = dat_tipic_wt, chains = n_chains, iter = n_iter, warmup = n_warmup, thin = n_thin)
fit_tipic_hf <- sampling(tstat_model, data = dat_tipic_hf, chains = n_chains, iter = n_iter, warmup = n_warmup, thin = n_thin)
fit_ticir_wt <- sampling(tstat_model, data = dat_ticir_wt, chains = n_chains, iter = n_iter, warmup = n_warmup, thin = n_thin)
fit_ticir_hf <- sampling(tstat_model, data = dat_ticir_hf, chains = n_chains, iter = n_iter, warmup = n_warmup, thin = n_thin)
fit_tpcpr_wt <- sampling(tstat_model, data = dat_tpcpr_wt, chains = n_chains, iter = n_iter, warmup = n_warmup, thin = n_thin)
fit_tpcpr_hf <- sampling(tstat_model, data = dat_tpcpr_hf, chains = n_chains, iter = n_iter, warmup = n_warmup, thin = n_thin)

# Save fits
save(ls(pattern='fit_'), file = 'C:/Users/Q/Dropbox/neon/code/tstats_stanfits.r')

# Diagnostic plots to see if models behave
# You can also plot pairs(fit, pars='beta') to see if any of the predictors are not giving additional information and should be dropped from the model.
traceplot(fit_tipic_wt, pars = 'beta')
stan_diag(fit_tipic_wt)
stan_rhat(fit_tipic_wt, pars = 'beta')

traceplot(fit_tipic_hf, pars = 'beta')
stan_diag(fit_tipic_hf)
stan_rhat(fit_tipic_hf, pars = 'beta')

traceplot(fit_ticir_wt, pars = 'beta')
stan_diag(fit_ticir_wt)
stan_rhat(fit_ticir_wt, pars = 'beta')

traceplot(fit_ticir_hf, pars = 'beta')
stan_diag(fit_ticir_hf)
stan_rhat(fit_ticir_hf, pars = 'beta')

traceplot(fit_tpcpr_wt, pars = 'beta')
stan_diag(fit_tpcpr_wt)
stan_rhat(fit_tpcpr_wt, pars = 'beta')

traceplot(fit_tpcpr_hf, pars = 'beta')
stan_diag(fit_tpcpr_hf)
stan_rhat(fit_tpcpr_hf, pars = 'beta')

# Extract summary info on parameters from model fit objects
summ_tipic_wt <- summary(fit_tipic_wt, probs=c(.025,.1,.5,.9,.975))
summ_tipic_hf <- summary(fit_tipic_hf, probs=c(.025,.1,.5,.9,.975))
summ_ticir_wt <- summary(fit_ticir_wt, probs=c(.025,.1,.5,.9,.975))
summ_ticir_hf <- summary(fit_ticir_hf, probs=c(.025,.1,.5,.9,.975))
summ_tpcpr_wt <- summary(fit_tpcpr_wt, probs=c(.025,.1,.5,.9,.975))
summ_tpcpr_hf <- summary(fit_tpcpr_hf, probs=c(.025,.1,.5,.9,.975))

# Display beta summary
summ_tipic_wt$summary[grep('beta', dimnames(summ_tipic_wt$summary)[[1]]), ]
summ_tipic_hf$summary[grep('beta', dimnames(summ_tipic_hf$summary)[[1]]), ]
summ_ticir_wt$summary[grep('beta', dimnames(summ_ticir_wt$summary)[[1]]), ]
summ_ticir_hf$summary[grep('beta', dimnames(summ_ticir_hf$summary)[[1]]), ]
summ_tpcpr_wt$summary[grep('beta', dimnames(summ_tpcpr_wt$summary)[[1]]), ]
summ_tpcpr_hf$summary[grep('beta', dimnames(summ_tpcpr_hf$summary)[[1]]), ]

# Reduced model
tstat_model_red1 <- stan_model('code/analysis/tstat_multiple_reg_reduced1.stan')
