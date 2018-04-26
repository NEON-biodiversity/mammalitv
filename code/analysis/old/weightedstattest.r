# Test weighted versus unweighted t-stats

source('code/analysis/Tstatsedited.r')

tsweighted <- QTstats(traits=traits_i, ind.plot = factor(mam_capture_sitemerge$plotID[i_i]), sp = factor(mam_capture_sitemerge$taxonID[i_i]), reg.pool=pool_i, nperm=5)

tsunweighted <- Tstats(traits=traits_i, ind.plot = factor(mam_capture_sitemerge$plotID[i_i]), sp = factor(mam_capture_sitemerge$taxonID[i_i]), reg.pool=pool_i, nperm=5)

cbind(tsweighted$Tstats$T_IP.IC[,'logweight'], tsunweighted$Tstats$T_IP.IC[,'logweight'])

source('~/GitHub/NEON/code/analysis/tstat2longform.r')

sesweighted <- tstat2longform_ses(tsweighted, bysite = FALSE)
sesunweighted <- tstat2longform_ses(tsunweighted, bysite = FALSE)

names(sesweighted)[4] <- "t_weighted"
names(sesunweighted)[4] <- "t_unweighted"

library(ggplot2)
plotdat <- left_join(sesweighted[,c(1,2,3,4,7)], sesunweighted[,c(1,2,3,4,7)]) %>%
  filter(trait == 'logweight')
ggplot(plotdat, aes(x = t_weighted, y = t_unweighted)) +
  facet_grid(~ stat) + geom_point()