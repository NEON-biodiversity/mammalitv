# Check species mean CV of body size by richness and temperature.
# Also check all community level CV

# Edit: instead of CV use quantiles (IQR). Does not assume normality. 17 May 2017.

rm(list=ls())
source('code/vis/newloadplotdat_nullswap.r')

qprobs <- c(0.05, 0.95)

commcvs <- mam_forrich %>% 
  transform(logweight = log10(weight)) %>%
  group_by(siteID) %>%
  summarize(community_cv = sd(logweight,na.rm=T)/mean(logweight,na.rm=T),
			community_iqr = diff(quantile(logweight, probs = qprobs, na.rm = T)))

spcvs <- mam_forrich %>% 
  transform(logweight = log10(weight)) %>%
  group_by(siteID, taxonID) %>%
  summarize(species_cv = sd(logweight,na.rm=T)/mean(logweight,na.rm=T),
            species_iqr = diff(quantile(logweight, probs = qprobs, na.rm = T)),
            species_n = sum(!is.na(logweight)))

spcvs_mean <- spcvs %>% ungroup %>%
  group_by(siteID) %>%
  summarize(byspecies_cv = weighted.mean(species_cv, species_n, na.rm=T),
			byspecies_iqr = weighted.mean(species_iqr, species_n, na.rm=T))

cvdat <- o2015 %>%
  filter(trait=='logweight') %>%
  dplyr::select(siteID, chao1, bio1, bio6, ostat_norm) %>%
  left_join(commcvs) %>% left_join(spcvs_mean)

library(ggplot2)

pdf('C:/Users/Q/google_drive/NEON_EAGER/Figures/revisionfigs/cvfigures.pdf',height=5,width=5)
ggplot(cvdat, aes(x=chao1, y=community_cv)) + geom_point() + theme_minimal() + labs(x='richness',y='communitywide cv')
ggplot(cvdat, aes(x=chao1, y=byspecies_cv)) + geom_point() + theme_minimal() + labs(x='richness',y='species average cv')
ggplot(cvdat, aes(x=bio6, y=community_cv)) + geom_point() + theme_minimal() + labs(x='min temp',y='communitywide cv')
ggplot(cvdat, aes(x=bio6, y=byspecies_cv)) + geom_point() + theme_minimal() + labs(x='min temp',y='species average cv')
dev.off()

# Added 05 May: improve plot themes, and include statistics.

# is the communitywide cv what is driving overlap?
ggplot(cvdat, aes(x=community_90, y=ostat_norm)) + geom_point() + theme_john + labs(y='overlap',x='communitywide cv') # yes.
ggplot(cvdat, aes(x=byspecies_90, y=ostat_norm)) + geom_point() + theme_john + labs(y='overlap',x='species average cv') # no.
ggplot(cvdat, aes(x=byspecies_90, y=qlogis(ostat_norm))) + geom_point() + theme_john + labs(y='overlap',x='species average cv') # no.

# "niche sorting" : is the overlap lower at sites where the niche sorting z-score is higher?
cvdat$niche_zscore <- subset(o2015sn, trait=='logweight')$ostat_norm_localnull_ses

ggplot(cvdat, aes(x=niche_zscore, y=ostat_norm)) + geom_point() + theme_john + labs(y='overlap',x='niche sorting score') # no.

# test correlations
lm_mechs <- lm(I(qlogis(ostat_norm)) ~ community_iqr * byspecies_iqr * niche_zscore, data=cvdat %>% filter(complete.cases(.)), na.action='na.pass')
library(MuMIn)
dredge(lm_mechs)
lm1 <- lm(I(qlogis(ostat_norm)) ~ community_iqr + byspecies_iqr, data=cvdat %>% filter(complete.cases(.)), na.action='na.pass')
summary(lm1)
lm_maineffects <- lm(I(qlogis(ostat_norm)) ~ community_iqr + byspecies_iqr + niche_zscore, data=cvdat %>% filter(complete.cases(.)), na.action='na.pass')
lm_maineffects_standard <- lm(I(qlogis(ostat_norm)) ~ scale(community_iqr) + scale(byspecies_iqr) + scale(niche_zscore), data=cvdat %>% filter(complete.cases(.)), na.action='na.pass')
summary(lm_maineffects_standard)
confint(lm_maineffects_standard)

lm_mechs <- lm(I(qlogis(ostat_norm)) ~ community_iqr * byspecies_iqr, data=cvdat, na.action='na.pass')
library(MuMIn)
dredge(lm_mechs)
lm1 <- lm(I(qlogis(ostat_norm)) ~ community_iqr + byspecies_iqr, data=cvdat, na.action='na.pass')
summary(lm1)

lm1std <- lm(I(qlogis(ostat_norm)) ~ scale(community_iqr) + scale(byspecies_iqr), data=cvdat, na.action='na.pass')
summary(lm1std)

# Clean up figures

reg1 <- lm(I(qlogis(ostat_norm)) ~ community_iqr, data = cvdat)
coef1 <- reg1$coefficients
reg2 <- lm(I(qlogis(ostat_norm)) ~ byspecies_iqr, data = cvdat)
coef2 <- reg2$coefficients

logitline <- function(b0, b1, x) plogis(b0 + b1 * x)

plot1 <- ggplot(cvdat, aes(x=community_iqr, y=ostat_norm)) + geom_point(size = 3) + 
  stat_function(geom='line', fun = logitline, args=list(b0 = coef1[1], b1 = coef1[2]), color = 'black', size = 0.8, n=9999) +
  geom_text(data=data.frame(community_iqr=Inf, ostat_norm=Inf, letter='a'), aes(label=letter), size=10, hjust=1.1, vjust=1.1) +
  theme_john + labs(y='Overlap',x='Community-wide body mass range') 
plot2 <- ggplot(cvdat, aes(x=byspecies_iqr, y=ostat_norm)) + geom_point(size = 3) + 
  stat_function(geom='line', fun = logitline, args=list(b0 = coef2[1], b1 = coef2[2]), color = 'black', size = 0.8, n=9999) +
  geom_text(data=data.frame(byspecies_iqr=Inf, ostat_norm=Inf, letter='b'), aes(label=letter), size=10, hjust=1.1, vjust=1.1) +
  theme_john + labs(y='Overlap',x='Species average body mass range')

library(gridExtra)
fp <- 'C:/Users/Q/google_drive/NEON_EAGER/Figures/revisionfigs'
png(file.path(fp, 'fig5mechanisms.png'), height=5, width=10, res=400, units='in')
grid.arrange(plot1, plot2, nrow=1)
dev.off()



# 24 July: better multiple-regression visualizations ----------------------

# https://stats.stackexchange.com/questions/125561/what-does-an-added-variable-plot-partial-regression-plot-explain-in-a-multiple

library(car)

avPlots(lm1std)
avPlots(lm_maineffects_standard)
avdat <- avPlots(lm1std)
avdat <- lapply(avdat, function(x) {x <- as.data.frame(x); names(x) <- c('x','y'); x})

logitline <- function(b0, b1, x) plogis(b0 + b1 * x)

avreg1 <- lm(y ~ x, data = avdat[[1]])
avcoef1 <- avreg1$coefficients
avreg2 <- lm(y ~ x, data = avdat[[2]])
avcoef2 <- avreg2$coefficients

avplot1 <- ggplot(avdat[[1]], aes(x = x, y = plogis(y))) + geom_point(size = 3) + 
  stat_function(geom='line', fun = logitline, args=list(b0 = avcoef1[1], b1 = avcoef1[2]), color = 'black', size = 0.8, n=9999) +
  geom_text(data=data.frame(x=Inf,y=Inf, letter='a'), aes(label=letter), size=10, hjust=1.1, vjust=0.9) +
  theme_john + labs(y='Overlap',x='Community-wide body mass range')

avplot2 <- ggplot(avdat[[2]], aes(x = x, y = plogis(y))) + geom_point(size = 3) + 
  stat_function(geom='line', fun = logitline, args=list(b0 = avcoef2[1], b1 = avcoef2[2]), color = 'black', size = 0.8, n=9999) +
  geom_text(data=data.frame(x=Inf,y=Inf, letter='b'), aes(label=letter), size=10, hjust=1.1, vjust=0.9) +
  theme_john + labs(y='Overlap',x='Species average body mass range')

library(gridExtra)
fp <- 'C:/Users/Q/google_drive/NEON_EAGER/Figures/msfigs2017jan/'
png(file.path(fp, 'fig5addedvariable_shortlabel.png'), height=5, width=10, res=400, units='in')
grid.arrange(avplot1, avplot2, nrow=1)
dev.off()

pdf(file.path(fp, 'fig5.pdf'), height=5, width=10)
grid.arrange(avplot1, avplot2, nrow=1)
dev.off()