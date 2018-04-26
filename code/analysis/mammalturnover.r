# Quantify temporal turnover between mammal surveys and relate to climate
# Author: QDR
# Project: NEON ITV
# Created: 18 July 2016
# Last modified: 21 July 2016 (expanded data added)

data_path <- '/mnt/research/neon'
# Calculate distance between consecutive mammal measurements or some other metric of beta-diversity


load(file.path(data_path, 'final_data/allorganismal2016Jun29.r')) # New organismal data 
load(file.path(data_path, 'final_data/allsiteplot2016Jun21.r')) # Site covariates

# Tabulate species by site and date

# Fix some of the poorly formatted dates in the table.
mam_capture <- subset(mam_capture, !is.na(siteID) & siteID!='')

library(lubridate)
#dateformats <- guess_formats(mam_capture$date, orders=c('mdY','Ymd'))
#sapply(head(mam_capture$date), function(x) ifelse(!is.na(ymd(x)), ymd(x), mdy(x)))

date1 <- ymd(mam_capture$date)
date1[is.na(date1)] <- mdy(mam_capture$date[is.na(date1)])
mam_capture$date <- date1

comm_array <- with(mam_capture, table(taxonID, date, siteID))

library(plyr)

# Sum up the individual nights into bouts (assign manually)

bouts <- list(c(1,1,1,2,2,2,3,3,3,4,4,4,4,4,4,5,5,5,5,5),
              c(1,2,2,2,3,3,3,3,3,3,4,4,4,4,5,5,5,5,6,6,6,7,7,7,7,8,8,8,9,9,9),
              c(1,1,1,2,3,4,4,4,5,5,5,6,6,6,7,7,7),
              c(1,1,1,2,2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4,5,5,5,5,5,5,6,6,7,7,7,8,8,8,8,8,8,9,9,9,9,9,9,10,10,10,10,10,10,11,11,11,11,11,11),
              c(1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,4,4,4,4),
              c(1,1,1),
              c(1,1,1,1,1,1,2,2,2,3,3,3,3,3,3,4,4,4,5,5,6,6,7,7,7,7,7,7,7),
              c(1,1,1,1,2,2,3,3,3,4,4,4,5,5,5,6,6,6,6,7,7,7,7,8,8,8,8),
              c(1,1,1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7,7,8,8,8,9,9,9,10,10,10,11,11,11,11,12),
              c(1,1,1,2,2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4,5,5,5,5,5,5,6,6,7,7,7,7,7),
              c(1,1,1,1,1,1),
              c(1,1,2,2,2,2,2,3,3,3,3,4,4,4,5,5,5),
              c(1,1,1,1,2,2,2,3,3,3,3,3,3,4,4,4,5,5,5,5,5,5,6,6,6,6,6,7,7,7,7),
              c(1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,3),
              c(1,1,1,1,1,1,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4,5,5,5,5,5,5,6,6,6,6,6,6))

comm_list <- dlply(mam_capture, .(siteID), function(x) table(x$date, x$taxonID) )

# Remove sites with only one bout.

comm_list <- comm_list[-c(6,11)]
bouts <- bouts[-c(6,11)]

# Calculate normalized relative abundance by bout for each site.

comm_bybout <- list()

for (i in 1:length(comm_list)) {
  x <- as.data.frame.matrix(comm_list[[i]])
  x$bout <- bouts[[i]]
  x$date <- as.Date(row.names(x))
  comm_bybout[[i]] <- ddply(x, .(bout), function(z) {
    zcomm <- z[, -c(ncol(z)-1, ncol(z))]
    abunds <- apply(zcomm, 2, sum)
    abunds <- abunds/sum(abunds)
    return(data.frame(date = z$date[1], t(abunds)))
  })
}

# Calculate the beta diversity at the different sites.

library(vegetarian)
# d(comm_bybout[[1]][,-(1:2)], lev='beta', q=1, boot=TRUE)

betadiv_all <- sapply(comm_bybout, function(x) as.numeric(d(x[,-(1:2)], lev='beta', q=1, boot=TRUE)))

betadiv_all <- as.data.frame(t(betadiv_all))
names(betadiv_all) <- c('D.Value','StdErr')
betadiv_all$siteID <- names(comm_list)
betadiv_all <- merge(betadiv_all, neonsitedata, sort=FALSE, all.x=TRUE, all.y=FALSE)

# Statistical tests, weighted by inverse of standard error
bio1lm <- lm(D.Value ~ bio1, weights = 1/StdErr, data = betadiv_all)
bio4lm <- lm(D.Value ~ bio4, weights = 1/StdErr, data = betadiv_all)
bio15lm <- lm(D.Value ~ bio15, weights = 1/StdErr, data = betadiv_all) # "Significant"

library(ggplot2)
ggplot(betadiv_all, aes(x=bio1, y=D.Value, ymin=D.Value-StdErr, ymax=D.Value+StdErr)) +
  geom_pointrange() +
  theme_minimal() +
  ggtitle('Within-year mammal turnover vs site temperature') +
  labs(x = 'Mean annual temperature (C)', y = 'Temporal beta-diversity')
ggsave('figs/png/mammalturnover_vs_mat.png', height=6, width=6, dpi=300)

ggplot(betadiv_all, aes(x=bio4, y=D.Value, ymin=D.Value-StdErr, ymax=D.Value+StdErr)) +
  geom_pointrange() +
  theme_minimal() +
  ggtitle('Within-year mammal turnover vs within-year temp variability') +
  labs(x = 'Temperature Seasonality (CV)', y = 'Temporal beta-diversity')
ggsave('figs/png/mammalturnover_vs_tempCV.png', height=6, width=6, dpi=300)


ggplot(betadiv_all, aes(x=bio15, y=D.Value)) +
  geom_pointrange(aes(ymin=D.Value-StdErr, ymax=D.Value+StdErr)) +
  theme_minimal() +
#  geom_line(aes(y = predict(bio15lm)), size=1.5) +
  ggtitle('Within-year mammal turnover vs within-year precip variability') +
 # geom_text(data=data.frame(bio15 = 40, D.Value = 1.65, lab = "list(R^2 == 0.43, p == 0.008)"), aes(label=lab), parse=TRUE) +
 # geom_text(data=data.frame(bio15 = 19, D.Value = 1.49, lab = "Outlier:\nOak Ridge"), aes(label=lab), color='red', hjust=.2) +
  labs(x = 'Precipitation Seasonality (CV)', y = 'Temporal beta-diversity')
ggsave('figs/png/mammalturnover_vs_precipCV.png', height=6, width=6, dpi=300)

# Another way of calculating would be to look at the distance/time between successive bouts.
# In other words, the rate of change at a given time of the community.

turnover_rates <- list()

for (i in 1:length(comm_bybout)) {
  z <- comm_bybout[[i]]
  turnover_rates[[i]] <- ddply(z, .(year(date)), function(x) {
    if (nrow(x) < 2) return(NA)
    days_i <- as.numeric(diff(as.Date(x$date)))
    rates_i <- c()
    for (j in 2:nrow(x)) {
      diff_j <- x[j, -(1:2)] - x[j-1, -(1:2)]
      euc_j <- sqrt(sum(diff_j^2))
      rates_i[j-1] <- euc_j/days_i[j-1]
    }
    mean(rates_i)
  })
}

for (i in 1:length(turnover_rates)) turnover_rates[[i]]$siteID <- names(comm_list)[i]
turnover_rates <- do.call('rbind', turnover_rates)
names(turnover_rates) <- c('year', 'variability', 'siteID')
rates2014 <- subset(turnover_rates, year==2014)

rates2014 <- merge(rates2014, neonsitedata, sort=FALSE, all.x=TRUE, all.y=FALSE)

library(ggplot2)
ggplot(rates2014, aes(x=bio1, y=variability)) +
  geom_point() +
  theme_minimal() +
  ggtitle('Within-year mammal turnover (method 2, 2014 only) vs site temperature') +
  labs(x = 'Mean annual temperature (C)', y = 'Temporal beta-diversity')
ggsave('figs/png/mammalturnover_method2_vs_mat.png', height=6, width=6, dpi=300)

ggplot(rates2014, aes(x=bio4, y=variability)) +
  geom_point() +
  theme_minimal() +
  ggtitle('Within-year mammal turnover vs within-year temp variability') +
  labs(x = 'Temperature Seasonality (CV)', y = 'Temporal beta-diversity')
ggsave('figs/png/mammalturnover_method2_vs_tempCV.png', height=6, width=6, dpi=300)




###############################

# beetle turnover

# Tabulate species by site and date

bet_pinning <- subset(bet_pinning, !is.na(siteID) & siteID!='')

# Correct dates that are not on the same bout in the bet_pinning data frame.

bet_pinning$collectDate[bet_pinning$siteID=='OSBS' & bet_pinning$collectDate=='2014-06-07'] <- '2014-06-09'
bet_pinning$collectDate[bet_pinning$siteID=='DSNY' & bet_pinning$collectDate=='2014-06-06'] <- '2014-06-11'
bet_pinning$collectDate[bet_pinning$siteID=='JERC' & bet_pinning$collectDate=='2014-09-10'] <- '2014-09-11'
bet_pinning$collectDate[bet_pinning$siteID=='ORNL' & bet_pinning$collectDate=='2014-06-17'] <- '2014-06-16'
bet_pinning$collectDate[bet_pinning$siteID=='SCBI' & bet_pinning$collectDate=='2014-06-11'] <- '2014-06-10'
bet_pinning$collectDate[bet_pinning$siteID=='SCBI' & bet_pinning$collectDate=='2014-08-21'] <- '2014-08-19'
bet_pinning$collectDate[bet_pinning$siteID=='UNDE' & bet_pinning$collectDate=='2014-08-12'] <- '2014-08-19'
bet_pinning$collectDate[bet_pinning$siteID=='WOOD' & bet_pinning$collectDate=='2014-08-17'] <- '2014-08-13'
bet_pinning$collectDate[bet_pinning$siteID=='WOOD' & bet_pinning$collectDate=='2014-08-23'] <- '2014-08-27'


beetlecomm_array <- with(bet_pinning, table(taxonID, collectDate, siteID))

library(plyr)

beetlecomm_list <- dlply(bet_pinning, .(siteID), function(x) table(x$collectDate, x$taxonID) )

# Remove the ones with no entries
beetlecomm_list <- beetlecomm_list[!names(beetlecomm_list) %in% c('ONAQ','TALL')]

# Calculate the beta diversity at the different sites.

library(vegetarian)
#d(comm_list[[1]], lev='beta', q=1, boot=TRUE)

beetlebetadiv_all <- sapply(beetlecomm_list, function(x) as.numeric(d(x, lev='beta', q=1, boot=TRUE)))

beetlebetadiv_all <- as.data.frame(t(beetlebetadiv_all))
names(beetlebetadiv_all) <- c('D.Value','StdErr')
beetlebetadiv_all$siteID <- names(beetlecomm_list)
beetlebetadiv_all <- merge(beetlebetadiv_all, neonsitedata, sort=FALSE, all.x=TRUE, all.y=FALSE)

# Statistical tests, weighted by inverse of standard error
bio1lm <- lm(D.Value ~ bio1, weights = 1/StdErr, data = beetlebetadiv_all) # "significant"
bio4lm <- lm(D.Value ~ bio4, weights = 1/StdErr, data = beetlebetadiv_all) # "significant"
bio15lm <- lm(D.Value ~ bio15, weights = 1/StdErr, data = beetlebetadiv_all)

ggplot(beetlebetadiv_all, aes(x=bio1, y=D.Value)) +
  geom_pointrange(aes(ymin=D.Value-StdErr, ymax=D.Value+StdErr)) +
  theme_minimal() +
  geom_line(aes(y=predict(bio1lm)), size=1.5) +
  geom_text(data=data.frame(bio1 = 12.5, D.Value = 3.5, lab = "list(R^2 == 0.50, p == 0.016)"), aes(label=lab), parse=TRUE) +
  qSubtitle('Temporal beetle beta-diversity vs site temperature', 'Weighted linear regression fit') +
  labs(x = 'Mean annual temperature (C)', y = 'Temporal beta-diversity')
ggsave('figs/png/beetleturnover_vs_mat.png', height=6, width=6, dpi=300)

ggplot(beetlebetadiv_all, aes(x=bio4, y=D.Value)) +
  geom_pointrange(aes(ymin=D.Value-StdErr, ymax=D.Value+StdErr)) +
  theme_minimal() +
  geom_line(aes(y=predict(bio4lm)), size=1.5) +
  geom_text(data=data.frame(bio4 = 875, D.Value = 3.5, lab = "list(R^2 == 0.49, p == 0.016)"), aes(label=lab), parse=TRUE) +
  qSubtitle('Within-year beetle turnover vs within-year temp variability', 'Weighted linear regression fit') +
  labs(x = 'Temperature Seasonality (CV)', y = 'Temporal beta-diversity')
ggsave('figs/png/beetleturnover_vs_tempCV.png', height=6, width=6, dpi=300)


ggplot(beetlebetadiv_all, aes(x=bio15, y=D.Value)) +
  geom_pointrange(aes(ymin=D.Value-StdErr, ymax=D.Value+StdErr)) +
  theme_minimal() +
  ggtitle('Within-year beetle turnover vs within-year precip variability') +
  labs(x = 'Precipitation Seasonality (CV)', y = 'Temporal beta-diversity')
ggsave('figs/png/beetleturnover_vs_precipCV.png', height=6, width=6, dpi=300)