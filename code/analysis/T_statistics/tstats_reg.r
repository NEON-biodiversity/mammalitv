# Basic regressions on T-statistics
# Author: QDR
# Project: NEON ITV
# Created: 19 Aug 2016
# Last modified: 24 Aug 2016

# Modified 24 Aug: fix incorrect thing on mammal grids, also add new site/plot covariates with the missing Alaska data

# Run regressions on the T-statistics, using the predictors specified in the grant proposal.
# For now, I am just using multiple regression. Eventually, we will use more sophisticated
# models if something interesting seems to be coming out of this.

# Mammal data will be used since it's the main one we have right now.


# Data loading/processing -------------------------------------------------

load('C:/Users/Q/Dropbox/neon/code/tstatsmam.r')
source('code/analysis/tstat2longform.r')

# Convert t-statistics outputs to long form for analysis.
tm_long14 <- tstat2longform_ses(tstats_mam_regpools14iucn)
tm_long15 <- tstat2longform_ses(tstats_mam_regpools15iucn)

# Load neon data, calculate species richness, load spatial stats, and merge everything with the T-statistics data frames.
load('allorganismal_latest.r')
load('allsiteplot_latest.r')


tm_clim14 <- merge(tm_long14, neonsitedata)
tm_clim15 <- merge(tm_long15, neonsitedata)

library(lubridate)
library(plyr)
library(reshape2)
#mam_capture <- transform(mam_capture, gridID = substr(trapCoordinate,1,1))

chao1mammal14 <- ddply(subset(mam_capture, year(date) == 2014), .(siteID), function(x) {
  xmat <- dcast(x, formula = plotID ~ taxonID)[,-1] # Get rid of date col, use number of rows as abundance
  S_obs <- ncol(xmat)
  f1 <- sum(apply(xmat, 2, sum) == 1)
  f2 <- sum(apply(xmat, 2, sum) == 2)
  return(data.frame(chao1 = S_obs + (f1 * (f1 - 1)) / (2 * (f2 + 1))))
})

tm_clim14 <- merge(tm_clim14, chao1mammal14)

chao1mammal15 <- ddply(subset(mam_capture, year(date) == 2015), .(siteID), function(x) {
  xmat <- dcast(x, formula = plotID ~ taxonID)[,-1] # Get rid of date col, use number of rows as abundance
  S_obs <- ncol(xmat)
  f1 <- sum(apply(xmat, 2, sum) == 1)
  f2 <- sum(apply(xmat, 2, sum) == 2)
  return(data.frame(chao1 = S_obs + (f1 * (f1 - 1)) / (2 * (f2 + 1))))
})

tm_clim15 <- merge(tm_clim15, chao1mammal15)

spatialstat <- read.csv('C:/Users/Q/Dropbox/neon/data/external_datasets/NEON_spatial_stats.csv')
neonplotdata <- cbind(neonplotdata, spatialstat)

library(dplyr)
spatialmeans <- neonplotdata %>% group_by(siteID) %>% 
  dplyr::summarize(ruggedness = mean(tri, na.rm=TRUE), roughness = mean(roughness, na.rm=TRUE)) %>% as.data.frame

tm_clim14 <- merge(tm_clim14, spatialmeans)
tm_clim15 <- merge(tm_clim15, spatialmeans)


# Multiple regression -----------------------------------------------------

# T_ip.ic
# Naive model (multiple regression with everything thrown in)
library(MuMIn)

# Include two-way interactions, use 2015 data only, and only use body weight and hind foot length because they have the most data.

# Weight
tipic_wtlmdat <- subset(tm_clim15, stat=='T_IP.IC' & !is.na(t) & !is.na(bio1) & trait=='weight', select=c(t, bio1, bio4, cv_bio1, chao1, ruggedness))

# Standardize the predictor variables
tipic_wtlmdat <- mutate_each(tipic_wtlmdat, funs((. - mean(.))/sd(.)), bio1:ruggedness)

# Estimate full linear model to 2nd order interactions
tipic_wtlm <- lm(t ~ .*., data = tipic_wtlmdat, na.action = 'na.pass')

# "Stupid" model selection to find best model
tipic_wtdredge <- dredge(tipic_wtlm)

head(tipic_wtdredge)

tipic_wtbest <- get.models(tipic_wtdredge, subset = delta < 5)

summary(tipic_wtbest[[1]])
hist(tipic_wtbest[[1]]$residuals) # Not bad but one large residual
# The only model that outperforms the null model is a richness-only model. As richness increases, Tipic increases. This is the opposite of prediction 1 from the grant.

# Hind Foot Length
tipic_hflmdat <- subset(tm_clim15, stat=='T_IP.IC' & !is.na(t) & !is.na(bio1) & trait=='hindfootLength', select=c(t, bio1, bio4, cv_bio1, chao1, ruggedness))

# Standardize the predictor variables
tipic_hflmdat <- mutate_each(tipic_hflmdat, funs((. - mean(.))/sd(.)), bio1:ruggedness)

# Estimate full linear model to 2nd order interactions
tipic_hflm <- lm(t ~ .*., data = tipic_hflmdat, na.action = 'na.pass')

# "Stupid" model selection to find best model
tipic_hfdredge <- dredge(tipic_hflm)

head(tipic_hfdredge)

# Null model is the best, but delta of ruggedness-only model is 1.
tipic_hfbest <- get.models(tipic_hfdredge, subset = delta < 5)

# Although this isn't significant, the coefficient is positive, which is what Hypothesis 3 predicts.
summary(tipic_hfbest[[2]])
hist(tipic_hfbest[[2]]$residuals)
plot(tipic_hfbest[[2]])

# T_ic.ir

# Weight
ticir_wtlmdat <- subset(tm_clim15, stat=='T_IC.IR' & !is.na(t) & !is.na(bio1) & trait=='weight', select=c(t, bio1, bio4, cv_bio1, chao1, ruggedness))

# Standardize the predictor variables
ticir_wtlmdat <- mutate_each(ticir_wtlmdat, funs((. - mean(.))/sd(.)), bio1:ruggedness)

# Estimate full linear model to 2nd order interactions
ticir_wtlm <- lm(t ~ .*., data = ticir_wtlmdat, na.action = 'na.pass')

# "Stupid" model selection to find best model
ticir_wtdredge <- dredge(ticir_wtlm)

head(ticir_wtdredge)

ticir_wtbest <- get.models(ticir_wtdredge, subset = delta < 5)

summary(ticir_wtbest[[1]])
hist(ticir_wtbest[[1]]$residuals)
plot(ticir_wtbest[[1]])
# This is a lot more interesting. There is a lot more variation explained here.

# Hind Foot Length
ticir_hflmdat <- subset(tm_clim15, stat=='T_IC.IR' & !is.na(t) & !is.na(bio1) & trait=='hindfootLength', select=c(t, bio1, bio4, cv_bio1, chao1, ruggedness))

# Standardize the predictor variables
ticir_hflmdat <- mutate_each(ticir_hflmdat, funs((. - mean(.))/sd(.)), bio1:ruggedness)

# Estimate full linear model to 2nd order interactions
ticir_hflm <- lm(t ~ .*., data = ticir_hflmdat, na.action = 'na.pass')

# "Stupid" model selection to find best model
ticir_hfdredge <- dredge(ticir_hflm)

head(ticir_hfdredge)

ticir_hfbest <- get.models(ticir_hfdredge, subset = delta < 5)

summary(ticir_hfbest[[1]])
hist(ticir_hfbest[[1]]$residuals)
plot(ticir_hfbest[[1]])
# Positive correlation with mean annual temperature and with species richness.





# Extra code below. Don't use.

# 2015 TIPIC
tipicdat15 <- tm_clim15 %>% filter(stat == 'T_IP.IC', !is.na(t)) %>% dplyr::select(t, trait, bio1, bio4, cv_bio1, chao1, ruggedness) 

# Standardize
tipicdat15 <- mutate_each(tipicdat15, funs((. - mean(.))/sd(.)), bio1:ruggedness)

tipiclm15 <- tipicdat15 %>% group_by(trait) %>% do(mod = lm(t ~ ., data = dplyr::select(., -trait), na.action='na.pass'))

do(tipiclm15, print(summary(.$mod))) # No terms are significant.

# 2015 TICIR
ticirdat15 <- tm_clim15 %>% filter(stat == 'T_IC.IR', !is.na(t) & !is.na(bio1)) %>% dplyr::select(t, trait, bio1, bio4, cv_bio1, chao1, ruggedness) 

# Standardize
ticirdat15 <- mutate_each(ticirdat15, funs((. - mean(., na.rm=T))/sd(., na.rm=T)), bio1:ruggedness)

ticirlm15 <- ticirdat15 %>% group_by(trait) %>% do(mod = lm(t ~ ., data = dplyr::select(., -trait), na.action='na.pass'))

do(ticirlm15, print(summary(.$mod))) # No terms are significant.
