# Format tstats object to be plotted with ggplot2 (extract summary stats from null distributions)
# Author: QDR
# Project: NEON ITV
# Created: 10 Aug 2016
# Last modified: 26 Aug 2016

# Modified 26 Aug: added option to return plotID instead of siteID
# Modified 18 Aug: added an additional function that computes standardized effect size

tstat2longform <- function(ts, bysite=TRUE) {
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
  if (!bysite) {
    names(out)[2] <- 'plotID'
    out$siteID <- substr(out$plotID, 1, 4)
  }
  row.names(out) <- 1:nrow(out)
  out
}

# Modified version of the tstat2longform function that computes standardized effect sizes

tstat2longform_ses <- function(ts, bysite=TRUE) {
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
  if (!bysite) {
    names(out)[2] <- 'plotID'
    out$siteID <- substr(out$plotID, 1, 4)
  }
  row.names(out) <- 1:nrow(out)
  out
}

overlap2longform_ses <- function(ts, bysite=TRUE) {
  ses_ovl <- ses_overlap(ts)
  ns <- dim(ts$Tstats$T_IP.IC_nm)
  ntrait <- ns[2]
  nsite <- ns[3]
  out <- list()
  for (trait in 1:ntrait) {
    for (site in 1:nsite) {
      overlap_i <- ses_ovl$ses[site, trait]
      ci_min_i <- ses_ovl$ses_lower[site, trait]
      ci_max_i <- ses_ovl$ses_upper[site, trait]
      ci_min_raw_i <- ses_ovl$raw_lower[site, trait]
      ci_max_raw_i <- ses_ovl$raw_upper[site, trait]
      out[[length(out)+1]] <- data.frame(trait = dimnames(ts$Tstats[[1]])[[2]][trait],
                                         siteID = dimnames(ts$Tstats[[1]])[[1]][site],
                                         overlap = ts$overlaps[site, trait],
                                         overlap_ses = overlap_i,
                                         ci_min = ci_min_i,
                                         ci_max = ci_max_i,
                                         ci_min_raw = ci_min_raw_i,
                                         ci_max_raw = ci_max_raw_i
      )
    }
  }
  
  out <- do.call('rbind', out)
  if (!bysite) {
    names(out)[2] <- 'plotID'
    out$siteID <- substr(out$plotID, 1, 4)
  }
  row.names(out) <- 1:nrow(out)
  out
}

