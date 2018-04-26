# SEM with productivity included
# QDR 28 March 2017

# Modified 05 April 2017: correct NPP (use log10-transformation instead of raw value)
# Modified 30 March 2017: add in some extra variables (try out different ones)

source('code/vis/newloadplotdat.r')
data_path <- '/mnt/research/neon'

modisdat <- read.csv(file.path(data_path, 'external_data/final_external_data/NEON_modisyearly_withcv.csv'), stringsAsFactors = FALSE)
modis_mean <- modisdat %>% group_by(siteID) %>% summarize_at(-(1:3), mean, na.rm=TRUE)

semdat <- o2015 %>%
  filter(trait == 'logweight') %>%
  left_join(modis_mean) %>%
  mutate(o = qlogis(ostat_norm),
         t = (bio1 - mean(bio1, na.rm=T))/sd(bio1, na.rm=T),
         r = (chao1 - mean(bio1, na.rm=T))/sd(chao1, na.rm=T),
         p = (NPP - mean(NPP, na.rm=T))/sd(NPP, na.rm=T)) %>%
  dplyr::select(siteID, o, t, r, p) 

library(GGally)
ggpairs(semdat[,-1])

full_mediation_model2 <- 
  ' # direct effect
r ~ e*t
# mediator
o ~ a*t + b*p
p ~ d*t
r ~ c*o + f*p
# indirect effect (a*c)
ac := a*c
dbc := d*b*c
df := d*f
# total effect
total := e + (a*c) + (d*b*c) + (d*f)
# intercepts
o ~ 1
r ~ 1
p ~ 1
'

library(blavaan)

mediation_bayes2 <- blavaan(model = full_mediation_model2, data = semdat,
                           auto.var = TRUE, auto.fix.first = TRUE, auto.cov.lv.x = TRUE,
                           jagcontrol = list(method = 'rjparallel'),
                           n.chains = 3)

parameterEstimates(mediation_bayes2)
inspect(mediation_bayes2, 'rsquare')

# Temperature and LAI are not at all related. Remove that path from the model.
reduced_mediation_model2 <- 
  ' # direct effect
r ~ e*t
# mediator
o ~ a*t + b*p
r ~ c*o + f*p
# indirect effect of temperature (a*c)
ac := a*c
# total effect of temperature
totaltemp := e + (a*c)
# indirect effect of productivity (b*c)
bc := b*c
# total effect of productivity
totalprod := f + (b*c)
# intercepts
o ~ 1
r ~ 1
'

set.seed(27609)

mediation_reduced <- blavaan(model = reduced_mediation_model2, data = semdat,
                            auto.var = TRUE, auto.fix.first = TRUE, auto.cov.lv.x = TRUE,
                            jagcontrol = list(method = 'rjparallel'),
                            n.chains = 3, burnin = 5000, sample = 20000)

parameterEstimates(mediation_reduced)
inspect(mediation_reduced, 'rsquare')

parameterEstimates(mediation_reduced, level = 0.90)

fitMeasures(mediation_bayes2)
fitMeasures(mediation_reduced)

# Get pseudo p-values out of the credible interval samples
# Concatenate the three chains first

allchains <- do.call('rbind', mediation_reduced@external$runjags$mcmc)

apply(allchains>0,2,sum)/nrow(allchains)
apply(allchains<0,2,sum)/nrow(allchains)


direct_temp <- allchains[,1]
direct_prod <- allchains[,5]

indirect_temp <- allchains[,2] * allchains[,4]
indirect_prod <- allchains[,3] * allchains[,4]

net_temp <- allchains[,1] + allchains[,2] * allchains[,4]
net_prod <- allchains[,5] + allchains[,3] * allchains[,4]

sum(direct_temp > 0)/length(direct_temp) # pseudo p-value 0.15
sum(direct_prod > 0)/length(direct_prod) # pseudo p-value 0.22
sum(indirect_temp < 0)/length(indirect_temp) # pseudo p-value 0.02
sum(indirect_prod > 0)/length(indirect_prod) # pseudo p-value 0.02
sum(net_temp < 0)/length(net_temp) # pseudo p-value 0.31
sum(net_prod > 0)/length(net_prod) # pseudo p-value 0.002


###############################
# Modification 30 March 2017: model selection to show that direct effects aren't important

nodirecteffects <- 
  ' # direct effects
#r ~ e*t
#r ~ f*p
# mediator
o ~ a*t + b*p
r ~ c*o
# indirect effect of temperature (a*c)
ac := a*c
# total effect of temperature
#totaltemp := e + (a*c)
# indirect effect of productivity (b*c)
bc := b*c
# total effect of productivity
#totalprod := f + (b*c)
# intercepts
o ~ 1
r ~ 1
'

directeffect_temp <-
  ' # direct effects
r ~ e*t
#r ~ f*p
# mediator
o ~ a*t + b*p
r ~ c*o
# indirect effect of temperature (a*c)
ac := a*c
# total effect of temperature
totaltemp := e + (a*c)
# indirect effect of productivity (b*c)
bc := b*c
# total effect of productivity
#totalprod := f + (b*c)
# intercepts
o ~ 1
r ~ 1
'

directeffect_prod <-
  ' # direct effects
#r ~ e*t
r ~ f*p
# mediator
o ~ a*t + b*p
r ~ c*o
# indirect effect of temperature (a*c)
ac := a*c
# total effect of temperature
#totaltemp := e + (a*c)
# indirect effect of productivity (b*c)
bc := b*c
# total effect of productivity
totalprod := f + (b*c)
# intercepts
o ~ 1
r ~ 1
'

set.seed(2865935)

nodirecteffects_fit <- blavaan(model = nodirecteffects, data = semdat,
                               auto.var = TRUE, auto.fix.first = TRUE, auto.cov.lv.x = TRUE,
                               jagcontrol = list(method = 'rjparallel'),
                               n.chains = 3, burnin = 5000, sample = 20000)
directeffect_temp_fit <- blavaan(model = directeffect_temp, data = semdat,
                                 auto.var = TRUE, auto.fix.first = TRUE, auto.cov.lv.x = TRUE,
                                 jagcontrol = list(method = 'rjparallel'),
                                 n.chains = 3, burnin = 5000, sample = 20000)
directeffect_prod_fit <- blavaan(model = directeffect_prod, data = semdat,
                                 auto.var = TRUE, auto.fix.first = TRUE, auto.cov.lv.x = TRUE,
                                 jagcontrol = list(method = 'rjparallel'),
                                 n.chains = 3, burnin = 5000, sample = 20000)

summary(nodirecteffects_fit)
parameterEstimates(nodirecteffects_fit)
inspect(nodirecteffects_fit, 'rsquare')
fitMeasures(nodirecteffects_fit)

summary(directeffect_temp_fit)
parameterEstimates(directeffect_temp_fit)
inspect(directeffect_temp_fit, 'rsquare')
fitMeasures(directeffect_temp_fit)

summary(directeffect_prod_fit)
parameterEstimates(directeffect_prod_fit)
inspect(directeffect_prod_fit, 'rsquare')
fitMeasures(directeffect_prod_fit)



#######################################

# 30 March: try different variables.
# temperature variables: bio1 (mat), bio4 (temp seasonality), bio6 (lowest monthly mean temp of the year)
# precipitation variables: bio12 (map), bio15 (precip seasonality), bio14 (lowest monthly precip total of the year)
# productivity variables: LAI, NDVI, NPP, cv_LAI, cv_NDVI

# z transform bio1, bio4, bio6, bio12, bio15, bio14, LAI, NDVI, NPP

fullsemdat <- o2015 %>%
  filter(trait == 'logweight') %>%
  left_join(modis_mean) %>%
  mutate(o = qlogis(ostat_norm),
         t = (bio1 - mean(bio1, na.rm=T))/sd(bio1, na.rm=T),
         r = (chao1 - mean(bio1, na.rm=T))/sd(chao1, na.rm=T),
         l = (LAI - mean(LAI, na.rm=T))/sd(LAI, na.rm=T),
         p = (decimalLongitude - mean(decimalLongitude, na.rm=T))/sd(decimalLongitude, na.rm=T),
         tvar = (bio4 - mean(bio4, na.rm=T))/sd(bio4, na.rm=T),
         tmin = (bio6 - mean(bio6, na.rm=T))/sd(bio6, na.rm=T),
         prec = (bio12 - mean(bio12, na.rm=T))/sd(bio12, na.rm=T),
         precvar = (bio15 - mean(bio15, na.rm=T))/sd(bio15, na.rm=T),
         precmin = (bio14 - mean(bio14, na.rm=T))/sd(bio14, na.rm=T),
         ndvi = (NDVI - mean(NDVI, na.rm=T))/sd(NDVI, na.rm=T)) %>%
  dplyr::select(siteID, o, t, r, p, l, tvar, tmin, prec, precvar, precmin, ndvi, cv_LAI, cv_NDVI) 

library(GGally)
ggpairs(fullsemdat[,-1])


library(blavaan)

bothdirecteffects <- 
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

tempdirecteffect <- 
  ' # direct effects
r ~ e*tmin
#r ~ f*p
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
#total_prod := f + (b*c)
# intercepts
o ~ 1
r ~ 1
'

proddirecteffect <- 
  ' # direct effects
#r ~ e*tmin
r ~ f*p
# mediator
o ~ a*tmin + b*p
r ~ c*o
# indirect effect of temperature (a*c)
indir_temp := a*c
# total effect of temperature
#total_temp := e + (a*c)
# indirect effect of productivity (b*c)
indir_prod := b*c
# total effect of productivity
total_prod := f + (b*c)
# intercepts
o ~ 1
r ~ 1
'

nodirecteffects <- 
  ' # direct effects
#r ~ e*tmin
#r ~ f*p
# mediator
o ~ a*tmin + b*p
r ~ c*o
# indirect effect of temperature (a*c)
indir_temp := a*c
# total effect of temperature
#total_temp := e + (a*c)
# indirect effect of productivity (b*c)
indir_prod := b*c
# total effect of productivity
#total_prod := f + (b*c)
# intercepts
o ~ 1
r ~ 1
'

set.seed(6885949)

bothdirect_fit <- blavaan(model = bothdirecteffects, data = fullsemdat,
                            auto.var = TRUE, auto.fix.first = TRUE, auto.cov.lv.x = TRUE,
                            jagcontrol = list(method = 'rjparallel'),
                            n.chains = 3, burnin = 5000, sample = 20000)

tempdirect_fit <- blavaan(model = tempdirecteffect, data = fullsemdat,
                          auto.var = TRUE, auto.fix.first = TRUE, auto.cov.lv.x = TRUE,
                          jagcontrol = list(method = 'rjparallel'),
                          n.chains = 3, burnin = 5000, sample = 20000)

proddirect_fit <- blavaan(model = proddirecteffect, data = fullsemdat,
                          auto.var = TRUE, auto.fix.first = TRUE, auto.cov.lv.x = TRUE,
                          jagcontrol = list(method = 'rjparallel'),
                          n.chains = 3, burnin = 5000, sample = 20000)

nodirect_fit <- blavaan(model = nodirecteffects, data = fullsemdat,
                          auto.var = TRUE, auto.fix.first = TRUE, auto.cov.lv.x = TRUE,
                          jagcontrol = list(method = 'rjparallel'),
                          n.chains = 3, burnin = 5000, sample = 20000)


parameterEstimates(bothdirect_fit)
parameterEstimates(tempdirect_fit)
parameterEstimates(proddirect_fit)
parameterEstimates(nodirect_fit)

inspect(bothdirect_fit, 'rsquare')
inspect(tempdirect_fit, 'rsquare')
inspect(proddirect_fit, 'rsquare')
inspect(nodirect_fit, 'rsquare')

fitMeasures(bothdirect_fit)
fitMeasures(tempdirect_fit)
fitMeasures(proddirect_fit)
fitMeasures(nodirect_fit)

allchains <- do.call('rbind', bothdirect_fit@external$runjags$mcmc)

apply(allchains,2,mean)

apply(allchains>0,2,sum)/nrow(allchains)
apply(allchains<0,2,sum)/nrow(allchains)


direct_temp <- allchains[,1]
direct_prod <- allchains[,2]

indirect_temp <- allchains[,3] * allchains[,5]
indirect_prod <- allchains[,4] * allchains[,5]

net_temp <- allchains[,1] + allchains[,3] * allchains[,5]
net_prod <- allchains[,2] + allchains[,4] * allchains[,5]

sum(direct_temp > 0)/length(direct_temp) # pseudo p-value 0.12
sum(direct_prod > 0)/length(direct_prod) # pseudo p-value 0.13
sum(indirect_temp < 0)/length(indirect_temp) # pseudo p-value 0.007
sum(indirect_prod > 0)/length(indirect_prod) # pseudo p-value 0.17
sum(net_temp < 0)/length(net_temp) # pseudo p-value 0.24
sum(net_prod > 0)/length(net_prod) # pseudo p-value 0.07


##############################
# 07 April: just run model with temperature and overlap. Should prefer indirect-only rather than direct + indirect.

tempdi <- 
' # direct effect
r ~ beta1*tmin
# mediator
o ~ beta2*tmin
r ~ beta3*o
# indirect effect of temperature (beta2*beta3)
indir_temp := beta2 * beta3
# total effect of temperature
total_temp := beta1 + beta2 * beta3
# intercepts
o ~ 1
r ~ 1
'

tempi <- 
  ' # direct effect
#r ~ beta1*tmin
# mediator
o ~ beta2*tmin
r ~ beta3*o
# indirect effect of temperature (beta2*beta3)
indir_temp := beta2 * beta3
# total effect of temperature
#total_temp := beta1 + beta2 * beta3
# intercepts
o ~ 1
r ~ 1
'

tempd <- 
' # direct effect
r ~ beta1*tmin
# mediator
#o ~ beta2*tmin
#r ~ beta3*o
# indirect effect of temperature (beta2*beta3)
#indir_temp := beta2 * beta3
# total effect of temperature
#total_temp := beta1 + beta2 * beta3
# intercepts
o ~ 1
r ~ 1
'

set.seed(4003846)

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

parameterEstimates(tempdi_fit)
parameterEstimates(tempi_fit)
parameterEstimates(tempd_fit)

inspect(tempdi_fit, 'rsquare')
inspect(tempi_fit, 'rsquare')
inspect(tempd_fit, 'rsquare')

fitMeasures(tempdi_fit)
fitMeasures(tempi_fit)
fitMeasures(tempd_fit)

# Pseudo p-values
allchains <- do.call('rbind', tempdi_fit@external$runjags$mcmc)

apply(allchains,2,mean)

apply(allchains>0,2,sum)/nrow(allchains)
apply(allchains<0,2,sum)/nrow(allchains)

direct_temp <- allchains[,1]
indirect_temp <- allchains[,2] * allchains[,3]
net_temp <- allchains[,1] + allchains[,2] * allchains[,3]

sum(direct_temp > 0)/length(direct_temp) # pseudo p-value 0.06
sum(indirect_temp < 0)/length(indirect_temp) # pseudo p-value 0.009
sum(net_temp < 0)/length(net_temp) # pseudo p-value 0.37
