# Regressions using bioclim pca1 as a variable
# Look at normals, interannual variability, and seasonality?

# load neon data
load('allorganismal2016Jun29.r') # Organismal data
load('allsiteplot2016Jun21.r') # Site covariates

# add bioclim PCA

climsites <- !is.na(neonsitedata$bio1)

bioclimpca <- prcomp(neonsitedata[climsites,6:24], scale=TRUE) # 19 bioclim variables

summary(bioclimpca)
bioclimpca$rotation

png('figs/png/bioclimpca.png', height=6, width=6, res=300, units ='in')
biplot(bioclimpca, xlab = 'PC1 ~ 50% of variation', ylab = 'PC2 ~ 30% of variation', main = 'NEON Bioclim PCA')
text(x = c(3, 3, -4), y = c(-4.5, 4, 3), labels = c('Increasing temp', 'Increasing precip', 'Increasing\ntemp variation'), col = 'blue')
dev.off()

bcpca_data <- matrix(NA, nrow=nrow(neonsitedata), ncol=3)
bcpca_data[climsites, ] <- bioclimpca$x[,1:3]
dimnames(bcpca_data)[[2]] <- c('PCA1','PCA2','PCA3')

neonsitedata <- cbind(neonsitedata, as.data.frame(bcpca_data))

# PC1 ~50% of variation. High value means warmer and wetter.
# PC2 adds another ~30%. High value indicates low diurnal temperature range, and low precipitation variability
# PC3 adds another ~10% so the total of the first 3 is 89%. Loaded on mean temperature of wettest quarter and precipitation of warmest quarter (identifies very warm and wet sites)

# Beetle community turnover

# Similar to the temperature trend
beetlepca1lm <- lm(D.Value ~ PCA1, weights=1/StdErr, data=beetlebetadiv_all)
summary(beetlepca1lm)

# Plant phenology things


# The budbreak ~ PCA1 seems to be more closely related than pure temperature.
phase1onset_pca1lm <- lmer(day_onset ~ PCA1 + (1|siteID) + (1|taxonID), weights = n_onset, data = subset(phenotable_sd, n_onset > 2 & phase==1))
r.squaredGLMM(phase1onset_pca1lm)
confint(phase1onset_pca1lm, parm='PCA1', method='boot')

# However the flowering date ~ PCA1 seems to be less related.
phase4onset_pca1lm <- lmer(day_onset ~ PCA1 + (1|siteID) + (1|taxonID), weights = n_onset, data = subset(phenotable_sd, n_onset > 2 & phase==4))
r.squaredGLMM(phase4onset_pca1lm)
confint(phase4onset_pca1lm, parm='PCA1', method='boot')



# Mammals: Bergmann's rule, and variability

####### Bergmann's rule: ########
# Recreate mam_adults
mam_adults <- merge(mam_adults[,1:11], neonsitedata, sort=FALSE, all.x=TRUE)

# Regressions
# Species, random intercepts

library(lme4)
adult_byspecies <- lmer(weight ~ PCA1 + (1|taxonID), data = mam_adults)
confint(adult_byspecies, parm='PCA1')
summary(adult_byspecies)

library(MuMIn)
r.squaredGLMM(adult_byspecies) # Very low rsq, but similar to temperature.

adult_sitespecies <- lmer(weight ~ PCA1 + taxonID + (1|siteID), data=mam_adults)
confint(adult_sitespecies, parm='PCA1')

r.squaredGLMM(adult_sitespecies)


# Mixed model with species and family
# I believe this is the correct nesting notation.
adult_taxa <- lmer(weight ~ bio1 + (1|family/taxonID), data = mam_adults)
confint(adult_taxa, parm='bio1')
r.squaredGLMM(adult_taxa)

