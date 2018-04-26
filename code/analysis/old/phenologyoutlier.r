# find places where some of the taxon ID are blank

names(phe_status)

# Each individualID should only have a single taxonID associated with it.

library(plyr)
library(lubridate)
taxdiag <- ddply(phe_status, .(individualID), summarize, ntax=length(unique(taxonID)), nmonth=length(unique(month(date))))

# Check sensitivity to outlier ('BART')

# There aren't enough observations any more to fit the model once the outlier is removed.

# By mean annual temperature
phase4onset_bio1lm <- lmer(day_onset ~ bio1 + (1|siteID), weights = n_onset, data = subset(phenotable_sd, n_onset > 2 & phase==4 & siteID!='BART'))
r.squaredGLMM(phase4onset_bio1lm)
confint(phase4onset_bio1lm, parm='bio1', method='boot')

# By interannual temperature variability
phase4onset_cvbio1lm <- lmer(day_onset ~ cv_bio1 + (1|siteID), weights = n_onset, data = subset(phenotable_sd, n_onset > 2 & phase==4 & siteID!='BART'))
r.squaredGLMM(phase4onset_cvbio1lm)
confint(phase4onset_cvbio1lm, parm='cv_bio1', method='boot')

# By within-year temperature seasonality
phase4onset_bio4lm <- lmer(day_onset ~ bio4 + (1|siteID), weights = n_onset, data = subset(phenotable_sd, n_onset > 2 & phase==4 & siteID!='BART'))
r.squaredGLMM(phase4onset_bio4lm)
confint(phase4onset_bio4lm, parm='bio4', method='boot')

