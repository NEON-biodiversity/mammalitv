# Customized T-statistics plots using ggplot2 (prettier than the default visualizations)
# Author: QDR
# Project: NEON ITV
# Created: 10 Aug 2016
# Last modified: 18 Aug 2016

# Modified 18 Aug: converted raw T-values to the standardized effect sizes relative to the null models.
# Modified 15 Aug: added spatial statistics plots, need to modify to split by year also.

# Load an example dataset (mammals)

library(cati)
library(ggplot2)
library(extrafont)
library(plyr)
library(reshape2)

load('C:/Users/Q/Dropbox/neon/code/tstatsmam.r')

# Format tstats object to be plotted with ggplot2 (extract summary stats from null distributions)
# Write function for the purpose

tstat2longform <- function(ts) {
  ns <- dim(ts$Tstats$T_IP.IC_nm)
  nperm <- ns[1]
  ntrait <- ns[2]
  nsite <- ns[3]
  stat_out <- list()
  for (stat in 1:3) {
    for (trait in 1:ntrait) {
      for (site in 1:nsite) {
        t_i <- ts$Tstats[[stat]][site, trait]
        t_null <- ts$Tstats[[stat+3]][, trait, site]
        t_null <- t_null[!is.na(t_null)]
        t_quant <- sum(t_null < t_i)/length(t_null)
        null_quant <- quantile(t_null, probs = c(0.025, 0.05, 0.1, 0.25, 0.50, 0.75, 0.9, 0.95, 0.975), na.rm=TRUE)
        stat_out[[length(stat_out) + 1]] <- data.frame(trait = dimnames(ts$Tstats[[1]])[[2]][trait],
                                                     siteID = dimnames(ts$Tstats[[1]])[[1]][site],
                                                     stat = names(ts$Tstats)[stat],
                                                     t = t_i,
                                                     quantile = t_quant,
                                                     null.025 = null_quant[1],
                                                     null.050 = null_quant[2],
                                                     null.100 = null_quant[3],
                                                     null.250 = null_quant[4],
                                                     null.500 = null_quant[5],
                                                     null.750 = null_quant[6],
                                                     null.900 = null_quant[7],
                                                     null.950 = null_quant[8],
                                                     null.975 = null_quant[9])
      }
    }
  }
  out <- do.call('rbind', stat_out)
  row.names(out) <- 1:nrow(out)
  out
}

# Modified version of the tstat2longform function that computes standardized effect sizes

tstat2longform_ses <- function(ts) {
  ses_list <- list(T_IP.IC = ses(ts$Tstats$T_IP.IC, ts$Tstats$T_IP.IC_nm),
                   T_IC.IR = ses(ts$Tstats$T_IC.IR, ts$Tstats$T_IC.IR_nm),
                   T_PC.PR = ses(ts$Tstats$T_PC.PR, ts$Tstats$T_PC.PR_nm))
  ns <- dim(ts$Tstats$T_IP.IC_nm)
  ntrait <- ns[2]
  nsite <- ns[3]
  stat_out <- list()
  for (stat in 1:3) {
    for (trait in 1:ntrait) {
      for (site in 1:nsite) {
        t_i <- ses_list[[stat]]$ses[site, trait]
        ci_min_i <- ses_list[[stat]]$ses.inf[site, trait]
        ci_max_i <- ses_list[[stat]]$ses.sup[site, trait]
        stat_out[[length(stat_out) + 1]] <- data.frame(trait = dimnames(ts$Tstats[[1]])[[2]][trait],
                                                       siteID = dimnames(ts$Tstats[[1]])[[1]][site],
                                                       stat = names(ts$Tstats)[stat],
                                                       t = t_i,
                                                       ci_min = ci_min_i,
                                                       ci_max = ci_max_i
                                                       )
      }
    }
  }
  out <- do.call('rbind', stat_out)
  row.names(out) <- 1:nrow(out)
  out
}


# Run function on example data
tm_long <- tstat2longform(tstats_mam_regpools)

# All-country regional species pool
tm_long_1pool <- tstat2longform(tstats_mam)

# Plot the t-statistics versus null distributions

ptstat <- ggplot(tm_long, aes(x=siteID)) + facet_grid(stat ~ trait, scales = 'free_y') +
  geom_segment(aes(y = null.025, yend = null.975, xend=siteID), alpha = 0.5, size = 1.5, color = 'skyblue') +
  geom_point(aes(y = t)) +
  labs(y = 'T-statistic') +
  theme_minimal() + theme(panel.border = element_rect(fill='transparent'))

ptstat + coord_flip()

# Log scale
ptstat + 
  geom_hline(yintercept = 1, linetype = 'dotted', color = 'midnightblue') +
  scale_y_log10() +
  coord_flip() + 
  labs(y = parse(text = 'paste(\"log\"[10], \" T-statistic\")'))
  

# Plot the t-statistics versus climate variables.

tm_clim <- merge(tm_long, neonsitedata)

source('code/bioclimnames.r')

pttemp <- ggplot(subset(tm_clim, trait %in% c('hindfootLength', 'weight', 'totalLength')), aes(x=bio1)) + 
  facet_grid(stat ~ trait, scales = 'free_y') +
  geom_segment(aes(y = null.025, yend = null.975, xend=bio1), alpha = 0.5, size = 1.5, color = 'skyblue') +
  geom_point(aes(y = t), size = 1.5) +
  labs(y = 'T-statistic', x = parse(text=bioclimnames[1])) +
  theme_minimal() + theme(panel.border = element_rect(fill='transparent'))

pttemp

# Log scale
source('~/qutil.r')

pttemp +
  geom_hline(yintercept = 1, linetype = 'dotted', color = 'midnightblue', size=0.9) +
  scale_y_log10() +
  labs(y = parse(text = 'paste(\"log\"[10], \" T-statistic\")')) +
  qSubtitle('T-statistics for NEON mammal traits', 'on log scale with 95% CI of null distributions')

ggsave('figs/png/tstats/mammaltstats_logscale_vsMAT.png', height=9, width=9, dpi=300)

# Plot the t-statistics versus richness.

# Mammals have abundance so we use the Chao1 estimator.

library(lubridate)
mam_capture <- transform(mam_capture, gridID = substr(trapCoordinate,1,1))


date1 <- ymd(mam_capture$date)
date1[is.na(date1)] <- mdy(mam_capture$date[is.na(date1)])
mam_capture$date <- date1

chao1mammal <- ddply(subset(mam_capture, year(date) == 2014), .(siteID), function(x) {
  xmat <- dcast(x, formula = gridID ~ taxonID)[,-1] # Get rid of date col, use number of rows as abundance
  S_obs <- ncol(xmat)
  f1 <- sum(apply(xmat, 2, sum) == 1)
  f2 <- sum(apply(xmat, 2, sum) == 2)
  return(data.frame(chao1 = S_obs + (f1 * (f1 - 1)) / (2 * (f2 + 1))))
})

tm_clim <- merge(tm_clim, chao1mammal)

ptchao <- ggplot(subset(tm_clim, trait %in% c('hindfootLength', 'weight', 'totalLength')), aes(x=chao1)) + 
  facet_grid(stat ~ trait, scales = 'free_y') +
  geom_segment(aes(y = null.025, yend = null.975, xend=chao1), alpha = 0.5, size = 1.5, color = 'skyblue') +
  geom_point(aes(y = t), size = 1.5) +
  labs(y = 'T-statistic', x = 'Chao1 Richness Estimator') +
  theme_minimal() + theme(panel.border = element_rect(fill='transparent'))

ptchao +
  geom_hline(yintercept = 1, linetype = 'dotted', color = 'midnightblue', size=0.9) +
  scale_y_log10() +
  labs(y = parse(text = 'paste(\"log\"[10], \" T-statistic\")')) +
  qSubtitle('T-statistics for NEON mammals vs. richness', 'on log scale with 95% CI of null distributions')

ggsave('figs/png/tstats/mammaltstats_logscale_vsrichness.png', height=9, width=9, dpi=300)



# Spatial Statistics ------------------------------------------------------

# Load spatial stats for plotting versus the Tstats.

spatialstat <- read.csv('C:/Users/Q/Dropbox/neon/data/external_datasets/NEON_spatial_stats.csv')

neonplotdata <- cbind(neonplotdata, spatialstat)

library(dplyr)

spatialmeans <- neonplotdata %>% group_by(siteID) %>% 
  summarize(ruggedness = mean(tri, na.rm=TRUE), roughness = mean(roughness, na.rm=TRUE))

neonsitedata <- left_join(neonsitedata, spatialmeans)
