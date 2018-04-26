# New structural equation models.
# Revisions of NEON MS
# QDR 28 March 2017

source('code/vis/newloadplotdat.r')

# logit transformation of overlap in SEM.

# 6 models.

# 1. no overlap effect, only direct temperature effect *** univariate regression
# 2. no temperature effect, only direct overlap effect *** univariate regression
# 3. only indirect temperature effect (via overlap) *** our "pet" theory
# 4. only direct temperature and overlap effects *** like a multiple regression
# 5. both direct and indirect temperature effects *** "full" model
# 6. null model, with no effect of either temperature or overlap

library(lavaan)

semdat <- o2015 %>% 
  filter(trait == 'logweight') %>%
  mutate(o = qlogis(ostat_norm),
         t = (bio1 - mean(bio1, na.rm=T))/sd(bio1, na.rm=T),
         r = (chao1 - mean(bio1, na.rm=T))/sd(chao1, na.rm=T)) %>%
  dplyr::select(siteID, o, t, r, mpd_z, pc1_productivityheterogeneity, pc2_topographyheterogeneity)

# Model specifications for the 6 models

model1 <-
  'r ~ t
   r ~ 1
   o ~ 1'

model2 <-
  'r ~ o
   r ~ 1'

model3 <-
  'r ~ o
   o ~ t
   r ~ 1
   o ~ 1'

model4 <-
  'r ~ o + t
   r ~ 1'

model5 <-
  'r ~ o + t
   o ~ t
   r ~ 1
   o ~ 1'

model6 <- 
  'r ~ 1
   o ~ 1'

# Fit SEMs
fit1 <- sem(model = model1, data = semdat)
fit2 <- sem(model = model2, data = semdat)
fit3 <- sem(model = model3, data = semdat)
fit4 <- sem(model = model4, data = semdat)
fit5 <- sem(model = model5, data = semdat)
fit6 <- sem(model = model6, data = semdat)

# Output information criteria for SEMs
fitMeasures(fit1)['bic']
fitMeasures(fit2)['bic']
fitMeasures(fit3)['bic']
fitMeasures(fit4)['bic']
fitMeasures(fit5)['bic']
fitMeasures(fit6)['bic']

# Output variation explained by SEMs
inspect(fit1, 'rsquare')
inspect(fit2, 'rsquare')
inspect(fit3, 'rsquare')
inspect(fit4, 'rsquare')
inspect(fit5, 'rsquare')
inspect(fit6, 'rsquare')

# Parameter estimates for SEMs
parameterEstimates(fit1)
parameterEstimates(fit2)
parameterEstimates(fit3)
parameterEstimates(fit4)
parameterEstimates(fit5)
parameterEstimates(fit6)

# Modification indices for SEMs
modificationIndices(fit1)
modificationIndices(fit2)
modificationIndices(fit3)
modificationIndices(fit4)
modificationIndices(fit5)
modificationIndices(fit6)

# Bayesian version
library(blavaan)

bfit3 <- blavaan(model = model3, data = semdat,
                 auto.var = TRUE, auto.fix.first = TRUE, auto.cov.lv.x = TRUE,
                 jagcontrol = list(method = 'rjparallel'),
                 n.chains = 3)

bfit4 <- blavaan(model = model4, data = semdat,
                 auto.var = TRUE, auto.fix.first = TRUE, auto.cov.lv.x = TRUE,
                 jagcontrol = list(method = 'rjparallel'),
                 n.chains = 3)

bfit5 <- blavaan(model = model5, data = semdat,
                 auto.var = TRUE, auto.fix.first = TRUE, auto.cov.lv.x = TRUE,
                 jagcontrol = list(method = 'rjparallel'),
                 n.chains = 3)


summary(bfit3)
summary(bfit4)
summary(bfit5)

fitMeasures(bfit3)
fitMeasures(bfit4)
fitMeasures(bfit5)

##########################

library(semPlot)

p3 <- semPlotModel(fit3)
p4 <- semPlotModel(fit4)
p5 <- semPlotModel(fit5)

semPaths(object = p3, what = 'std', intercepts = FALSE, thresholds = FALSE, sizeMan = 8, edge.label.cex = 1.5, label.cex = 2.5)
semPaths(object = p4, what = 'std', intercepts = FALSE, thresholds = FALSE, sizeMan = 8, edge.label.cex = 1.5, label.cex = 2.5)
semPaths(object = p5, what = 'std', intercepts = FALSE, thresholds = FALSE, sizeMan = 8, edge.label.cex = 1.5, label.cex = 2.5)



##########################

# Mediation model
# Source: http://lavaan.ugent.be/tutorial/mediation.html

# exogenous variable (X) is temperature, z-transformed 
# mediator (M) is overlap, logit-transformed
# endogenous variable (Y) is Chao1 richness, z-transformed

mediation_dat <- o2015 %>% 
  filter(trait == 'logweight') %>%
  mutate(M = qlogis(ostat_norm),
         X = (bio1 - mean(bio1, na.rm=T))/sd(bio1, na.rm=T),
         Y = (chao1 - mean(bio1, na.rm=T))/sd(chao1, na.rm=T)) %>%
  dplyr::select(X, M, Y)

full_mediation_model <- 
  ' # direct effect
    Y ~ c*X
    # mediator
    M ~ a*X
    Y ~ b*M
    # indirect effect (a*b)
    ab := a*b
    # total effect
    total := c + (a*b)
    # intercepts
    M ~ 1
    Y ~ 1
'

mediation_fit <- sem(full_mediation_model, data = mediation_dat)

summary(mediation_fit)

mediation_bayes <- blavaan(model = full_mediation_model, data = mediation_dat,
                           auto.var = TRUE, auto.fix.first = TRUE, auto.cov.lv.x = TRUE,
                           jagcontrol = list(method = 'rjparallel'),
                           n.chains = 3)

summary(mediation_bayes)
