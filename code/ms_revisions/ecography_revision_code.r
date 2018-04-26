# Additional R code needed to address reviewers' comments in Ecography
# QDR 07 Nov 2017 (NEON ITV)

# Edited 25 Apr 2018: Load files from hpcc, not dropbox
data_path <- '/mnt/research/neon'

# Rerun SEMs without cper, niwo, and moab ---------------------------------

# CPER, NIWO, and MOAB were identified by the iNEXT analysis as being poorly sampled so that we do not have a good idea of their true richness.
# Reviewer 1 asked to rerun the SEMs without them.

load(file.path(data_path, 'MS1_RodentOverlap/R_data/workspace_07nov.RData'))

library(dplyr)
library(blavaan)

sites_exclude_sem <- c('CPER', 'NIWO', 'MOAB')

fullsemdat <- ostats %>%
  mutate(o = qlogis(ostat_norm),
         r = (chao1 - mean(chao1, na.rm=T))/sd(chao1, na.rm=T),
         p = (log10(NPP) - mean(log10(NPP), na.rm=T))/sd(log10(NPP), na.rm=T),
         tmin = (bio6 - mean(bio6, na.rm=T))/sd(bio6, na.rm=T)) %>%
  select(siteID, o, r, p, tmin) 

exclsemdat <- filter(fullsemdat, !siteID %in% sites_exclude_sem)

# Model specifications for SEMs

# This is the full model with indirect (overlap-mediated) effects and direct effects for both temperature (bio6) and productivity (NPP)
tempdi_proddi <- 
  ' # direct effects
r ~ e*tmin
r ~ f*p
# mediator
o ~ a*tmin + b*p
r ~ c*o
# indirect effect of temperature (a*c)
indir_temp := a*c
# total effect of temperature
total_temp := e + (a*c)
# indirect effect of productivity (b*c)
indir_prod := b*c
# total effect of productivity
total_prod := f + (b*c)
# intercepts
o ~ 1
r ~ 1
'

# Indirect and direct temperature effects, with no productivity effects
tempdi <- 
  ' # direct effects
r ~ e*tmin
# mediator
o ~ a*tmin
r ~ c*o
# indirect effect of temperature (a*c)
indir_temp := a*c
# total effect of temperature
total_temp := e + (a*c)
# intercepts
o ~ 1
r ~ 1
'

# Indirect temperature effect, with no productivity effects
tempi <- 
  ' 
# mediator
o ~ a*tmin
r ~ c*o
# indirect effect of temperature (a*c)
indir_temp := a*c
# intercepts
o ~ 1
r ~ 1
'

# Direct temperature effect, with no productivity effects
tempd <- 
  ' # direct effects
r ~ e*tmin
# mediator
r ~ c*o
# intercepts
r ~ 1
'

# Direct and indirect temperature effects, only indirect productivity effect
tempdi_prodi <- 
  ' # direct effects
r ~ e*tmin
# mediator
o ~ a*tmin + b*p
r ~ c*o
# indirect effect of temperature (a*c)
indir_temp := a*c
# total effect of temperature
total_temp := e + (a*c)
# indirect effect of productivity (b*c)
indir_prod := b*c
# intercepts
o ~ 1
r ~ 1
'

# Indirect temperature effect, direct and indirect productivity effects
tempi_proddi <- 
  ' # direct effects
r ~ f*p
# mediator
o ~ a*tmin + b*p
r ~ c*o
# indirect effect of temperature (a*c)
indir_temp := a*c
# indirect effect of productivity (b*c)
indir_prod := b*c
# total effect of productivity
total_prod := f + (b*c)
# intercepts
o ~ 1
r ~ 1
'

# Indirect temperature effect, indirect productivity effect
tempi_prodi <- 
  ' 
# mediator
o ~ a*tmin + b*p
r ~ c*o
# indirect effect of temperature (a*c)
indir_temp := a*c
# indirect effect of productivity (b*c)
indir_prod := b*c
# intercepts
o ~ 1
r ~ 1
'

# Run SEM model fits

set.seed(4003846)

tempdi_proddi_fit <- blavaan(model = tempdi_proddi, data = fullsemdat,
                             auto.var = TRUE, auto.fix.first = TRUE, auto.cov.lv.x = TRUE,
                             jagcontrol = list(method = 'rjparallel'),
                             n.chains = 3, burnin = 5000, sample = 20000)
tempdi_fit <- blavaan(model = tempdi, data = fullsemdat,
                      auto.var = TRUE, auto.fix.first = TRUE, auto.cov.lv.x = TRUE,
                      jagcontrol = list(method = 'rjparallel'),
                      n.chains = 3, burnin = 5000, sample = 20000)
tempi_fit <- blavaan(model = tempi, data = fullsemdat,
                     auto.var = TRUE, auto.fix.first = TRUE, auto.cov.lv.x = TRUE,
                     jagcontrol = list(method = 'rjparallel'),
                     n.chains = 3, burnin = 5000, sample = 20000)
tempd_fit <- blavaan(model = tempd, data = fullsemdat,
                     auto.var = TRUE, auto.fix.first = TRUE, auto.cov.lv.x = TRUE,
                     jagcontrol = list(method = 'rjparallel'),
                     n.chains = 3, burnin = 5000, sample = 20000)
tempdi_prodi_fit <- blavaan(model = tempdi_prodi, data = fullsemdat,
                            auto.var = TRUE, auto.fix.first = TRUE, auto.cov.lv.x = TRUE,
                            jagcontrol = list(method = 'rjparallel'),
                            n.chains = 3, burnin = 5000, sample = 20000)
tempi_proddi_fit <- blavaan(model = tempi_proddi, data = fullsemdat,
                            auto.var = TRUE, auto.fix.first = TRUE, auto.cov.lv.x = TRUE,
                            jagcontrol = list(method = 'rjparallel'),
                            n.chains = 3, burnin = 5000, sample = 20000)
tempi_prodi_fit <- blavaan(model = tempi_prodi, data = fullsemdat,
                           auto.var = TRUE, auto.fix.first = TRUE, auto.cov.lv.x = TRUE,
                           jagcontrol = list(method = 'rjparallel'),
                           n.chains = 3, burnin = 5000, sample = 20000)

# Fit same models with the three sites removed.
set.seed(6880535)

extempdi_proddi_fit <- blavaan(model = tempdi_proddi, data = exclsemdat,
                             auto.var = TRUE, auto.fix.first = TRUE, auto.cov.lv.x = TRUE,
                             jagcontrol = list(method = 'rjparallel'),
                             n.chains = 3, burnin = 5000, sample = 20000)
extempdi_fit <- blavaan(model = tempdi, data = exclsemdat,
                      auto.var = TRUE, auto.fix.first = TRUE, auto.cov.lv.x = TRUE,
                      jagcontrol = list(method = 'rjparallel'),
                      n.chains = 3, burnin = 5000, sample = 20000)
extempi_fit <- blavaan(model = tempi, data = exclsemdat,
                     auto.var = TRUE, auto.fix.first = TRUE, auto.cov.lv.x = TRUE,
                     jagcontrol = list(method = 'rjparallel'),
                     n.chains = 3, burnin = 5000, sample = 20000)
extempd_fit <- blavaan(model = tempd, data = exclsemdat,
                     auto.var = TRUE, auto.fix.first = TRUE, auto.cov.lv.x = TRUE,
                     jagcontrol = list(method = 'rjparallel'),
                     n.chains = 3, burnin = 5000, sample = 20000)
extempdi_prodi_fit <- blavaan(model = tempdi_prodi, data = exclsemdat,
                            auto.var = TRUE, auto.fix.first = TRUE, auto.cov.lv.x = TRUE,
                            jagcontrol = list(method = 'rjparallel'),
                            n.chains = 3, burnin = 5000, sample = 20000)
extempi_proddi_fit <- blavaan(model = tempi_proddi, data = exclsemdat,
                            auto.var = TRUE, auto.fix.first = TRUE, auto.cov.lv.x = TRUE,
                            jagcontrol = list(method = 'rjparallel'),
                            n.chains = 3, burnin = 5000, sample = 20000)
extempi_prodi_fit <- blavaan(model = tempi_prodi, data = exclsemdat,
                           auto.var = TRUE, auto.fix.first = TRUE, auto.cov.lv.x = TRUE,
                           jagcontrol = list(method = 'rjparallel'),
                           n.chains = 3, burnin = 5000, sample = 20000)

# 8. Extract fit statistics and compare SEMs ------------------------------

all_sem_stats <- function(fit) {
  cat('Parameter estimates\n-------------------\n')
  print(parameterEstimates(fit))
  cat('\nR-square\n--------\n')
  print(inspect(fit, 'rsquare'))
  cat('\nFit measures\n------------\n')
  print(fitMeasures(fit))
  cat('\nPseudo-p-values\n---------------\n')
  
  allchains <- do.call('rbind', fit@external$runjags$mcmc)
  print(apply(allchains>0,2,sum)/nrow(allchains))
  print(apply(allchains<0,2,sum)/nrow(allchains))
  
}

all_sem_stats(tempdi_proddi_fit)
all_sem_stats(tempdi_fit)
all_sem_stats(tempi_fit)
all_sem_stats(tempd_fit)
all_sem_stats(tempdi_prodi_fit)
all_sem_stats(tempi_proddi_fit)
all_sem_stats(tempi_prodi_fit)

# Comparison of BICs from each model.
fitnames <- c('tempdi_proddi_fit', 'tempdi_fit', 'tempi_fit', 'tempd_fit', 'tempdi_prodi_fit', 'tempi_proddi_fit', 'tempi_prodi_fit')

thebics <- sapply(fitnames, function(x) {
  c(bic_full = fitMeasures(get(x))['bic'],
    bic_excl = fitMeasures(get(paste0('ex',x)))['bic'])
})

ther2s <- sapply(fitnames, function(x) {
  c(r2_full = inspect(get(x), 'rsquare')[1],
    r2_excl = inspect(get(paste0('ex',x)), 'rsquare')[1])
})

# 9. Calculate community- and species-level quantile ranges ---------------

qprobs <- c(0.05, 0.95)

commiqrs <- final_mammal_data %>% 
  transform(logweight = log10(weight)) %>%
  group_by(siteID) %>%
  summarize(community_iqr = diff(quantile(logweight, probs = qprobs, na.rm = T)))

spiqrs <- final_mammal_data %>% 
  transform(logweight = log10(weight)) %>%
  group_by(siteID, taxonID) %>%
  summarize(species_iqr = diff(quantile(logweight, probs = qprobs, na.rm = T)),
            species_n = sum(!is.na(logweight)))

spiqrs_mean <- spiqrs %>% ungroup %>%
  group_by(siteID) %>%
  summarize(byspecies_iqr = weighted.mean(species_iqr, species_n, na.rm=T))

iqrdat <- ostats %>%
  dplyr::select(siteID, chao1, bio1, bio6, ostat_norm) %>%
  left_join(commiqrs) %>% left_join(spiqrs_mean)

iqrdat$niche_zscore <- ostats$ostat_norm_localnull_ses

# 10. Run linear model for three drivers of overlap -----------------------

library(MuMIn)

lm_maineffects_standard <- lm(I(qlogis(ostat_norm)) ~ scale(community_iqr) + scale(byspecies_iqr) + scale(niche_zscore), data=iqrdat %>% filter(complete.cases(.)), na.action='na.pass')

dredge(lm_maineffects_standard, rank = 'BIC')
summary(lm_maineffects_standard)
confint(lm_maineffects_standard)


# Check plant richness as a predictor -------------------------------------

plantrichness <- read.csv(file.path(data_path, 'external_data/final_external_data/NEON_background_plant_richness.csv'), stringsAsFactors = FALSE)

with(left_join(fullsemdat, plantrichness), plot(p, richness_Ellis))
with(left_join(fullsemdat, plantrichness), plot(p, richness_NEON))

lmrichellis <- with(left_join(fullsemdat, plantrichness), lm(richness_Ellis ~ tmin))
lmrichneon <- with(left_join(fullsemdat, plantrichness), lm(richness_NEON ~ tmin))

lmrichellis <- with(left_join(fullsemdat, plantrichness), lm(r ~ richness_Ellis))
lmrichneon <- with(left_join(fullsemdat, plantrichness), lm(r ~ richness_NEON))

fullsemdat <- ostats %>%
  left_join(plantrichness) %>%
  mutate(o = qlogis(ostat_norm),
         r = (chao1 - mean(chao1, na.rm=T))/sd(chao1, na.rm=T),
         p = (log10(NPP) - mean(log10(NPP), na.rm=T))/sd(log10(NPP), na.rm=T),
         tmin = (bio6 - mean(bio6, na.rm=T))/sd(bio6, na.rm=T),
         rplant = (richness_Ellis - mean(richness_Ellis, na.rm=T))/sd(richness_Ellis, na.rm=T)) %>%
  select(siteID, o, r, p, tmin, rplant) 

# What are pairwise correlations?
cor(fullsemdat[,-1])

# Try to fit some SEMs with plant richness in as a term (in addition to, or instead of, productivity)
# Try temp -> plant richness -> overlap -> mammal richness
# Also try temp -> overlap -> mammal richness with an independent effect of plant richness on overlap?
fourlinkchain <- 
' 
# mediator
rplant ~ a*tmin
o ~ b*rplant
r ~ c*o
# indirect effect of plant richness (a*b*c)
indir_rplant := a*b*c
# indirect effect of temperature (b*c)
indir_temp := b*c
# intercepts
o ~ 1
r ~ 1
rplant ~ 1
'

fourlink_fit <- blavaan(model = fourlinkchain, data = fullsemdat,
                     auto.var = TRUE, auto.fix.first = TRUE, auto.cov.lv.x = TRUE,
                     jagcontrol = list(method = 'rjparallel'),
                     n.chains = 3, burnin = 5000, sample = 20000)

all_sem_stats(fourlink_fit)

indep_rplant <-
' 
# mediator
o ~ a*rplant + b*tmin
r ~ c*o
# indirect effect of plant richness on mammal richness (a*c)
indir_rplant := a*c
# indirect effect of temperature on mammal richness (b*c)
indir_temp := b*c
# intercepts
o ~ 1
r ~ 1
'

indep_rplant_fit <- blavaan(model = indep_rplant, data = fullsemdat,
                        auto.var = TRUE, auto.fix.first = TRUE, auto.cov.lv.x = TRUE,
                        jagcontrol = list(method = 'rjparallel'),
                        n.chains = 3, burnin = 5000, sample = 20000)

all_sem_stats(indep_rplant_fit)

# Add columns to supplements ----------------------------------------------

source('code/vis/newloadplotdat_nullswap.r')

# Appendix 1: add latitude, longitude, vegetation type (nlcd), and rodent richness.

# Everything but NLCD veg. type.
ostats %>% 
  dplyr::select(siteName, siteID, decimalLatitude, decimalLongitude, elevation, bio1, bio12, ruggedness, Observed, chao1)

# Load NLCD veg type
nlcdpoints <- read.csv(file.path(data_path, 'external_data/final_external_data/NEON_nlcdpoints.csv'))

lctables <- neonplotdata %>%
  cbind(nlcdpoints) %>%
  filter(subtype == 'mammalGrid') %>%
  group_by(siteID) %>%
  do(lct = table(.$landcover2011))

most_common_nlcd <- sapply(lctables$lct, function(x) ifelse(length(x) > 0, as.numeric(names(sort(x, decreasing=TRUE))[1]), NA))



nlcdlegend <- c('open_water'=11, 'perennial_ice_snow'=12, 'developed_open_space'=21, 'developed_low_intensity'=22, 'barren_land'=31, 'deciduous_forest'=41, 'evergreen_forest'=42, 'mixed_forest'=43, 'dwarf_scrub'=51, 'shrub_scrub'=52, 'grassland_herbaceous'=71, 'pasture_hay'=81, 'cultivated_crops'=82, 'woody_wetlands'=90, 'emergent_herbaceous_wetlands'=95)

site_supp_table <- ostats %>%
  left_join(data.frame(siteID = lctables$siteID, nlcd = names(nlcdlegend)[match(most_common_nlcd, nlcdlegend)])) %>%
  dplyr::select(siteName, siteID, decimalLatitude, decimalLongitude, elevation, bio1, bio12, ruggedness, Observed, chao1, nlcd)

write.csv(site_supp_table, file = 'C:/Users/Q/google_drive/NEON_EAGER/Manuscript1_RodentOverlap/revision/ecography_resub/site_supp_table.csv', row.names = FALSE)

# Appendix 2: add mean and standard deviation

mammal_supp_table <-
  mam_forrich %>% filter(Pineda_Main_food != 'Insectivore', !siteID %in% c('DSNY','DELA','HEAL'), year ==2015, !grepl('sp\\.', scientificName)) %>%
  group_by(scientificName) %>%
  summarize(n_indiv = length(unique(individualandtag)),
            n_sites = length(unique(siteID)),
            mass_mean = mean(weight, na.rm = T),
            mass_sd = sd(weight, na.rm = T),
            mass_q1 = quantile(weight, prob=0.25, na.rm=T),
            mass_q3 = quantile(weight, prob=0.75, na.rm=T))

write.csv(mammal_supp_table, file = 'C:/Users/Q/google_drive/NEON_EAGER/Manuscript1_RodentOverlap/revision/ecography_resub/mammal_supp_table.csv', row.names = FALSE)


# Edit color scheme of fig 2 map ------------------------------------------

rm(list=ls())
source('code/vis/newloadplotdat.r')

# Panel a: Map

datamap <- ggplot(subset(neonsitedata, decimalLatitude>20 & decimalLatitude<50), aes(x=decimalLongitude, y=decimalLatitude)) +
  borders('state', fill='transparent') +
  coord_map() + 
  scale_x_continuous(limits = c(-125,-66.2), expand=c(0,0), breaks = c(-120, -105, -90, -75), labels = c('120° W', '105° W', '90° W', '75° W')) +
  scale_y_continuous(limits=c(24.9,50), expand=c(0,0), breaks = c(25, 35, 45), labels = c('25° N', '35° N', '45° N')) +
  theme(panel.border = element_rect(color='black', fill=NA), panel.background=element_blank(), panel.grid=element_blank(), axis.title = element_blank())

mamcs <- scale_fill_gradientn(name = 'Overlap', colours = RColorBrewer::brewer.pal(9, 'Oranges'), breaks=c(0,.25,.5,.75,1))

us_map_ovl <- datamap + geom_point(size=5, shape=21, aes(fill = ostat_norm), data=subset(o2015, decimalLatitude>20 & decimalLatitude<50 & trait == 'logweight')) +
  mamcs +
  theme(legend.position=c(0.01,0.01), legend.justification=c(0,0), legend.direction='horizontal')

ggsave('C:/Users/Q/google_drive/NEON_EAGER/Manuscript1_RodentOverlap/revision/ecography_resub/fig2map_new.png', us_map_ovl, height=4.15, width=9*5/6, dpi=400)
