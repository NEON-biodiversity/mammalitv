# Tests of Bergmann's Rule for mammals
# Author: QDR
# Project: NEON ITV
# Created: 18 July 2016
# Last modified: 22 July 2016

# Modified 20 Apr 2018: Remove reference to local files

# 22 July: add families
# 21 July: add new data

# Load data
load('allorganismal2016Jun29.r') # New org. data
load('allsiteplot2016Jun21.r') # Site covariates

# Add families
# Load full taxonomic categories
mammalTax <- read.csv('species_lists/mammalSpeciesList.csv', stringsAsFactors = FALSE)
mam_capture <- merge(mam_capture, mammalTax[,c('taxonID', 'family')], sort=FALSE, all.x=TRUE, all.y=FALSE)

# Clean data
#mam_capture_subset <- subset(mam_capture, !siteID %in% c('DSNY'))
mam_capture_subset <- mam_capture # Don't remove Disney any more.
library(plyr)
mam_noID <- mam_capture_subset[mam_capture_subset$individualID=='', c('siteID','taxonID','family','individualID','hindfootLength','earLength','tailLength','totalLength','weight','sex','lifeStage')]
mam_byindiv <- ddply(subset(mam_capture_subset, individualID != ''),
                     .(siteID, taxonID, family, individualID), 
                     colwise(median, .(hindfootLength, earLength, tailLength, totalLength, weight), na.rm=TRUE))
# Get sex and life stage for each individual
mam_byindiv_class <- ddply(subset(mam_capture_subset, individualID != ''),
                          .(siteID, taxonID, family, individualID),
                          function(x) data.frame(sex = x$sex[1],
                                                 lifeStage = names(sort(table(x$lifeStage),decreasing=TRUE))[1]))

mam_byindiv <- cbind(mam_byindiv, mam_byindiv_class[,c('sex','lifeStage')])                          
mam_byindiv <-rbind(mam_byindiv, mam_noID)

#mam_capture_merge <- merge(mam_byindiv, neonplotdata[,-c(1:6,9:10)], sort=FALSE, all.x=TRUE)
mam_capture_sitemerge <- merge(mam_byindiv, neonsitedata[,-c(3:5)], sort=FALSE, all.x=TRUE)

# Test Bergmann's Rule

# Naive model
br_naive <- lm(weight ~ bio1, data = mam_capture_sitemerge)
# We see the opposite. Significantly positive slope.

# Model within species. Random intercept assumed for each taxon
library(lme4)
br_byspecies <- lmer(weight ~ bio1 + (1|taxonID), data = mam_capture_sitemerge)
confint(br_byspecies)
# We still see a positive slope (Significant at 0.05)

# Remove juveniles (reduces from ~2000 to ~1200)
mam_adults <- subset(mam_capture_sitemerge, lifeStage=='adult')
adult_naive <- lm(weight ~ bio1, data = mam_adults)
adult_byspecies <- lmer(weight ~ bio1 + (1|taxonID), data = mam_adults)
# This does not qualitatively change the pattern. It does tighten the confidence interval a bit.

adult_byfamily <- lmer(weight ~ bio1 + (1|family), data = mam_adults)
#adult_nestedtaxa <- lmer(weight ~ bio1 + (taxonID|family), data = mam_adults) # not working well, too many parameters.
adult_taxa <- lmer(weight ~ bio1 + (1|taxonID) + (1|family), data = mam_adults)

# Control for the effect of sex as well, using only adults
# Model sex as fixed effect
bysexandspecies <- lmer(weight ~ bio1 + sex + (1|taxonID), data = mam_adults)
# If this is correct, it also doesn't change the pattern.

# Site as random effect

adult_sitespecies <- lmer(weight ~ bio1 + (1|taxonID) + (1|siteID), data=mam_adults)
adult_sitespecies <- lmer(weight ~ bio1 + taxonID + (1|siteID), data=mam_adults)
adult_site <- lmer(weight ~ bio1 + (1|siteID), data=mam_adults)
adult_sitefamily <- lmer(weight ~ bio1 + family + (1|siteID), data=mam_adults)

# Plot Bergmann's Rule

library(ggplot2)
ggplot(mam_adults, aes(x = bio1, y = weight)) +
  geom_point(aes(color = taxonID)) +
  stat_smooth(aes(group = taxonID), method = 'lm', se = FALSE) +
  theme_minimal()

# Plot the results of the lm

newdat <- with(subset(mam_adults, !is.na(weight)), expand.grid(taxonID = unique(taxonID), bio1 = c(min(bio1), max(bio1))))
newdat2 <- with(subset(mam_adults, !is.na(weight)), expand.grid(siteID = unique(siteID), bio1 = c(min(bio1), max(bio1))))

source('~/qutil.r')

ggplot(subset(mam_adults, !is.na(weight)), aes(x = bio1, y = weight)) +
  geom_point(aes(color = taxonID)) +
  geom_line(aes(y = predict(adult_byspecies), group = taxonID)) +
  theme_minimal() +
  qSubtitle('Opposite of Bergmann\'s Rule', 'Linear mixed model fit with random intercepts') +
  labs(x = 'Mean annual temperature (C)', y = 'Body weight (g)')
ggsave('figs/png/notbergmannsrule.png', height=6, width=6, dpi=300)


ggplot(subset(mam_adults, !is.na(weight)), aes(x=bio1, y=weight)) +
  geom_point(aes(color=taxonID)) +
  geom_line(aes(y= predict(adult_sitespecies), group=taxonID)) +
  theme_minimal() +
  qSubtitle('Opposite of Bergmann\'s Rule', 'Species fixed effect, site random effect') +
  labs(x = 'Mean annual temperature (C)', y = 'Body weight (g)')
ggsave('figs/png/notbergmannsrule_siteandspecies.png', height=6, width=6, dpi=300)


ggplot(subset(mam_adults, !is.na(weight)), aes(x=bio1, y=weight)) +
  geom_point(aes(color=family)) +
  geom_line(aes(y= predict(adult_sitefamily), group=family)) +
  theme_minimal() +
  qSubtitle('Opposite of Bergmann\'s Rule', 'Taxonomic family fixed effect, site random effect') +
  labs(x = 'Mean annual temperature (C)', y = 'Body weight (g)')
ggsave('figs/png/notbergmannsrule_siteandfamily.png', height=6, width=6, dpi=300)


# Bergmann's rule tested with species means at sites only, ignoring intraspecific variation.

library(plyr)
mam_spmeans <- ddply(mam_adults, .(siteID, bio1, taxonID), function(x) with(x, data.frame(weight = mean(weight, na.rm=T), n = sum(!is.na(weight)))))
                    
mam_spmeans$taxonID <- factor(mam_spmeans$taxonID)
spmean_lm <- lm(weight ~ bio1, weights=n, data=mam_spmeans) # Similar coefficient but wider confidence interval
spmean_lmer <- lmer(weight ~ bio1 + (1|taxonID), data=mam_spmeans) # Again, similar coefficient but wider confidence interval.

# Look at slopes of individual species.
# Limit those with >= 10 individuals

commonspp <- names(table(mam_adults$taxonID)[table(mam_adults$taxonID) >= 10])
mam_common <- subset(mam_adults, taxonID %in% commonspp)

library(ggplot2)
ggplot(mam_common, aes(x = bio1, y = weight)) + 
  stat_smooth(method = 'lm') +
  geom_point(aes(color = siteID)) +
  facet_wrap(~ taxonID) +
  theme_minimal()
  

# Allen's Rule: "corollary" to Bergmann's rule, that individuals tend to have longer appendages in warmer environments.
ar_naive <- lm(hindfootLength ~ bio1, data=mam_adults) # Significantly positive slope but very low rsq
ar_byspecies <- lmer(hindfootLength ~ bio1 + (1|taxonID), data = mam_adults) # After giving species random intercepts, becomes very significantly negative. Strange.
ar_bysexandspecies <- lmer(hindfootLength ~ bio1 + sex + (1|taxonID), data = mam_adults)

# Also test "residual condition" of mammals

residcond <- lm(log(weight) ~ log(totalLength), data=mam_adults)$residuals
# Only 596 have both measured but this should be enough to look at patterns.

#residcond_footlength <- lm(log(weight) ~ log(hindfootLength), data=mam_adults)$residuals
# Do not use the foot length as it is not a good indicator of body condition.

mam_adults$residualCondition <- residcond[row.names(mam_adults)]
rc_naive <- lm(residualCondition ~ bio1, data=mam_adults) # Significantly less in warmer environments
rc_byspecies <- lmer(residualCondition ~ bio1 + (1|taxonID), data = mam_adults) # Appears to hold after giving species random intercepts as well
rc_bysexandspecies <- lmer(residualCondition ~ bio1 + sex + (1|taxonID), data = mam_adults) # Morphology seems different between sexes (not interesting)
rc_sitespecies <- lmer(residualCondition ~ bio1 + taxonID + (1|siteID), data = mam_adults)
