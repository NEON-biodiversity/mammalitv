# Corrected analysis and visualization of Bergmann's Rule plot
# Author: QDR
# Project: NEON ITV
# Created: 03 August 2016
# Last modified: 10 August 2016

# Modified 09 Aug, 10 Aug: added multipanel plots by species and family, and new analysis

# newdat <- expand.grid(family = unique(mam_adults$family), 
#                       taxonID = unique(mam_adults$taxonID)[1],
#                       bio1 = c(min(mam_adults$bio1), max(mam_adults$bio1)))
# 
# newdat <- ddply(subset(mam_adults, !is.na(weight)), .(family), function(x) with(x, expand.grid(taxonID = unique(x$taxonID), bio1 = c(min(x$bio1), max(x$bio1)))))
# 
# newdat <- ddply(subset(mam_adults, !is.na(weight)), .(family), function(x) data.frame(bio1 = c(min(x$bio1), max(x$bio1))))
# newdat <- ddply(subset(mam_adults, !is.na(weight)), .(family, taxonID), function(x) data.frame(bio1 = c(min(x$bio1), max(x$bio1))))

# Create data frame for new data prediction.

library(dplyr)
newdat <- mam_adults %>% filter(!is.na(weight) & siteID!='HEAL' & !is.na(family)) %>%
  group_by(family, taxonID) %>%
  do(data.frame(bio1 = c(min(.$bio1), max(.$bio1)),
                siteID = c(.$siteID[which.min(.$bio1)], .$siteID[which.max(.$bio1)])
                ))

newdat2 <- mam_adults %>% filter(!is.na(weight) & siteID!='HEAL' & !is.na(family)) %>%
  group_by(family) %>%
  do(data.frame(bio1 = c(min(.$bio1), max(.$bio1)),
                siteID = c(.$siteID[which.min(.$bio1)], .$siteID[which.max(.$bio1)])
  ))


adult_taxa <- lmer(weight ~ bio1 + (1|family/taxonID) + (1|siteID), data = subset(mam_adults, !is.na(weight) & siteID!='HEAL'))
confint(adult_taxa, parm='bio1', method='boot')
r.squaredGLMM(adult_taxa)

library(ggplot2)
source('~/qutil.r')
source('code/bioclimnames.r')

ggplot(subset(mam_adults, !is.na(weight))) +
  geom_point(aes(color=family, y=weight, x=bio1)) +
  geom_line(aes(y= predict(adult_taxa, newdata=newdat, type='response'), x=bio1, group=taxonID, color=family), data=newdat) +
  theme_minimal() +
  qSubtitle('Opposite of Bergmann\'s Rule', 'Random intercepts for sites and species nested within family') +
  labs(x = parse(text=bioclimnames[1]), y = 'Body weight (g)')
ggsave('figs/png/notbergmannsrule_correctslopes.png', height=5, width=5, dpi=300)

#################

# Just family, no taxonID
newdat <- ddply(subset(mam_adults, !is.na(weight)), .(family), function(x) data.frame(bio1 = c(min(x$bio1), max(x$bio1)), siteID = c(x$siteID[which.min(x$bio1)], x$siteID[which.max(x$bio1)])))

adult_justfamily <- lmer(weight ~ bio1 + (1|family) + (1|siteID), data = mam_adults)

ggplot(subset(mam_adults, !is.na(weight))) +
  geom_point(aes(color=family, y=weight, x=bio1)) +
  geom_line(aes(y= predict(adult_justfamily, newdata=newdat2), x=bio1, group=family, color=family), data=newdat2) +
  theme_minimal() +
  qSubtitle('Opposite of Bergmann\'s Rule', 'Random intercepts for sites and families') +
  labs(x = parse(text=bioclimnames[1]), y = 'Body weight (g)')

# Include both line types in a single figure.
ggplot(subset(mam_adults, !is.na(weight))) +
  geom_point(aes(color=family, y=weight, x=bio1), alpha=0.6) +
  geom_line(aes(y= predict(adult_taxa, newdata=newdat, type='response'), size='Species', x=bio1, group=taxonID, color=family), data=newdat, alpha=0.4) +
  geom_line(aes(y= predict(adult_taxa, newdata=mutate(newdat2, taxonID=NA), re.form=~(1|family)), size='Family', x=bio1, group=family, color=family), data=newdat2) +
  scale_size_manual(name='Predictions', values=c('Species'=0.5,'Family'=2)) +
  theme_minimal() +
  qSubtitle('Opposite of Bergmann\'s Rule', 'Random intercepts for sites and families') +
  labs(x = parse(text=bioclimnames[1]), y = 'Body weight (g)')

# Species-specific analysis -----------------------------------------------

# Subset mam_adults by species with at least 3 sites

library(plyr)

sptable <- with(mam_adults, table(taxonID, siteID))
nsites <- apply(sptable > 0, 1, sum)
nsites[nsites > 2]

mam_multisite <- subset(mam_adults, taxonID %in% names(nsites[nsites > 2]))

mam_multisite_cv <- ddply(mam_multisite, .(siteID, taxonID, family), colwise(raster::cv, .cols = .(hindfootLength,earLength,tailLength,totalLength,weight), na.rm=T))
mam_multisite_n <- ddply(mam_multisite, .(siteID, taxonID, family), colwise(function(x) sum(!is.na(x)), .cols = .(hindfootLength, earLength,tailLength,totalLength,weight)))
names(mam_multisite_n) <- paste(names(mam_multisite_n), 'n', sep='_')

mam_multisite_cv <- cbind(mam_multisite_cv, mam_multisite_n[,4:8])
mam_multisite_cv <- merge(mam_multisite_cv, neonsitedata)

# Visualization
ggplot(mam_multisite_cv, aes(x=bio1, y=weight, size=weight_n, color=taxonID)) + geom_point() + theme_minimal()

# Faceted plot
ggplot(mam_multisite_cv, aes(x=bio1, y=weight, size=weight_n)) +
  facet_wrap(~ taxonID) +
  labs(x = parse(text=bioclimnames[1]), y='Body weight coefficient of variation', size='No. individuals') +
  geom_point() + theme_bw()
ggsave('figs/png/mammalcvplots/weightCV_vs_mat_byspecies.png', height=8*2/3, width=10*2/3, dpi=300)

# Add a line showing the temperature range it covers for each species

mamranges <- ddply(mam_multisite_cv, .(taxonID), summarize, mintemp=min(bio1), maxtemp=max(bio1))
ssetmam <- mam_multisite_cv$taxonID[!is.na(mam_multisite_cv$weight)]
ggplot(subset(mam_multisite_cv, !is.na(weight)), aes(x=bio1, y=weight, size=weight_n)) +
  facet_wrap(~ taxonID, drop = TRUE) +
  labs(x = parse(text=bioclimnames[1]), y='Body weight coefficient of variation', size='No. individuals') +
  geom_point() + 
  geom_segment(aes(x = mintemp, xend = maxtemp, y = 0, yend = 0, color = maxtemp-mintemp), data = subset(mamranges, taxonID %in% ssetmam), size = 1.5) +
  scale_color_gradient(name = 'Temperature\nrange', low = 'goldenrod', high = 'midnightblue') +
  #geom_line(aes(y = 0), color = 'red', size=1.5) +
  theme_bw()
ggsave('figs/png/mammalcvplots/weightCV_vs_mat_byspecieswithrange.png', height=8*2/3, width=10*2/3, dpi=300)

ggplot(mam_multisite_cv, aes(x=bio1, y=weight, size=weight_n)) +
  facet_wrap(~ family) +
  labs(x = parse(text=bioclimnames[1]), y='Body weight coefficient of variation', size='No. individuals') +
  geom_point() + theme_bw()
ggsave('figs/png/mammalcvplots/weightCV_vs_mat_byfamily.png', height=8*2/3, width=10*2/3, dpi=300)


# Faceted plot by interannual temp variation
ggplot(mam_multisite_cv, aes(x=cv_bio1, y=weight, size=weight_n)) +
  facet_wrap(~ taxonID) +
  labs(x = 'Interannual temperature coefficient of variation', y='Body weight coefficient of variation', size='No. individuals') +
  geom_point() + theme_bw() + theme(axis.text.x = element_text(size=6))
ggsave('figs/png/mammalcvplots/weightCV_vs_tempcv_byspecies.png', height=8*2/3, width=10*2/3, dpi=300)

# Add color for niche width to this plot.
ggplot(merge(subset(mam_multisite_cv, !is.na(weight)), mamranges), aes(x=cv_bio1, y=weight, size=weight_n, color=maxtemp-mintemp)) +
  facet_wrap(~ taxonID) +
  labs(x = 'Interannual temperature coefficient of variation', y='Body weight coefficient of variation', size='No. individuals') +
  scale_color_gradient(name = 'Temperature\nrange', low = 'goldenrod', high = 'midnightblue') +
  geom_point() + theme_bw() + theme(axis.text.x = element_text(size=6))
ggsave('figs/png/mammalcvplots/weightCV_vs_tempcv_byspeciescolored.png', height=8*2/3, width=10*2/3, dpi=300)

ggplot(mam_multisite_cv, aes(x=cv_bio1, y=weight, size=weight_n)) +
  facet_wrap(~ family) +
  labs(x = 'Interannual temperature coefficient of variation', y='Body weight coefficient of variation', size='No. individuals') +
  geom_point() + theme_bw()
ggsave('figs/png/mammalcvplots/weightCV_vs_tempcv_byfamily.png', height=8*2/3, width=10*2/3, dpi=300)


# Hypothesis testing
library(lme4)
library(MuMIn)
wtcv_bio1lm <- lmer(weight ~ bio1 + (1|taxonID) + (1|siteID), weights = weight_n, data = mam_multisite_cv)
wtcv_bio1lmfamily <- lmer(weight ~ bio1 + (1|family) + (1|siteID), weights = weight_n, data = mam_multisite_cv)
wtcv_bio1cvlm <- lmer(weight ~ cv_bio1 + (1|taxonID) + (1|siteID), weights = weight_n, data = mam_multisite_cv)
wtcv_bio1cvlmfamily <- lmer(weight ~ cv_bio1 + (1|family) + (1|siteID), weights = weight_n, data = mam_multisite_cv)

summary(wtcv_bio1lm)
confint(wtcv_bio1lm, parm='bio1', method='boot', nsim=999)
r.squaredGLMM(wtcv_bio1lm)

summary(wtcv_bio1lmfamily)
confint(wtcv_bio1lmfamily, parm='bio1', method='boot', nsim=999)
r.squaredGLMM(wtcv_bio1lmfamily)

summary(wtcv_bio1cvlm)
confint(wtcv_bio1cvlm, parm='cv_bio1', method='boot', nsim=999)
r.squaredGLMM(wtcv_bio1cvlm)

summary(wtcv_bio1cvlmfamily)
confint(wtcv_bio1cvlmfamily, parm='cv_bio1', method='boot', nsim=999)
r.squaredGLMM(wtcv_bio1cvlmfamily)


# Species as fixed effect

wtcv_bio1_spfixed <- lmer(weight ~ bio1 + taxonID + (1|siteID), weights=weight_n, data=mam_multisite_cv)
summary(wtcv_bio1_spfixed)
confint(wtcv_bio1_spfixed, method = 'boot', nsim = 999)
r.squaredGLMM(wtcv_bio1_spfixed)

wtcv_bio1cv_spfixed <- lmer(weight ~ cv_bio1 + taxonID + (1|siteID), weights=weight_n, data=mam_multisite_cv)
summary(wtcv_bio1cv_spfixed)
confint(wtcv_bio1cv_spfixed, method = 'boot', nsim = 999)
r.squaredGLMM(wtcv_bio1cv_spfixed)