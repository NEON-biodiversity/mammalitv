# Include longitude in the SEM to account for biogeographical effect
# Rough proxy for proximity to the highly diverse sierra madre corridor in Mexico

data_path <- '/mnt/research/neon'
source('code/vis/newloadplotdat.r')

modisdat <- read.csv(file.path(data_path, 'external_data/final_external_data/NEON_modisyearly_withcv.csv'), stringsAsFactors = FALSE)
modis_mean <- modisdat %>% group_by(siteID) %>% summarize_at(-(1:3), mean, na.rm=TRUE)

# z transform bio1, bio4, bio6, bio12, bio15, bio14, LAI, NDVI, NPP

fullsemdat <- o2015 %>%
  filter(trait == 'logweight') %>%
  left_join(modis_mean) %>%
  mutate(o = qlogis(ostat_norm),
         t = (bio1 - mean(bio1, na.rm=T))/sd(bio1, na.rm=T),
         r = (chao1 - mean(bio1, na.rm=T))/sd(chao1, na.rm=T),
         l = (LAI - mean(LAI, na.rm=T))/sd(LAI, na.rm=T),
         p = (NPP - mean(NPP, na.rm=T))/sd(NPP, na.rm=T),
         tvar = (bio4 - mean(bio4, na.rm=T))/sd(bio4, na.rm=T),
         tmin = (bio6 - mean(bio6, na.rm=T))/sd(bio6, na.rm=T),
         prec = (bio12 - mean(bio12, na.rm=T))/sd(bio12, na.rm=T),
         precvar = (bio15 - mean(bio15, na.rm=T))/sd(bio15, na.rm=T),
         precmin = (bio14 - mean(bio14, na.rm=T))/sd(bio14, na.rm=T),
         ndvi = (NDVI - mean(NDVI, na.rm=T))/sd(NDVI, na.rm=T),
         long = (decimalLongitude - mean(decimalLongitude, na.rm=T))/sd(decimalLongitude, na.rm=T)) %>%
  dplyr::select(siteID, o, t, r, p, l, tvar, tmin, prec, precvar, precmin, ndvi, cv_LAI, cv_NDVI, long) 


library(blavaan)

directeffects_long <- 
  ' # direct effects
r ~ beta1 * tmin
r ~ beta2 * p
r ~ beta3 * long
# mediator
o ~ beta5 * tmin + beta6 * p + beta7 * long
r ~ beta4 * o
# indirect effect of temperature (beta4beta5)
indir_temp := beta5 * beta4
# total effect of temperature
total_temp := beta1 + (beta5 * beta4)
# indirect effect of productivity (beta6beta4)
indir_prod := beta6 * beta4
# total effect of productivity
total_prod := beta2 + (beta6 * beta4)
# indirect effect of longitude
indir_long := beta7 * beta4
# total effect of longitude
total_long := beta3 + (beta7 * beta4)
# intercepts
o ~ 1
r ~ 1
'

set.seed(8675309)

bothdirect_fit <- blavaan(model = directeffects_long, data = fullsemdat,
                          auto.var = TRUE, auto.fix.first = TRUE, auto.cov.lv.x = TRUE,
                          jagcontrol = list(method = 'rjparallel'),
                          n.chains = 3, burnin = 5000, sample = 20000)

parameterEstimates(bothdirect_fit)
inspect(bothdirect_fit, 'rsquare')
fitMeasures(bothdirect_fit)
