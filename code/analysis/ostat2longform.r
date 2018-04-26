# O-stat output to longform data frame.
# Author: QDR
# Project: NEON ITV
# Created: 3 Oct 2016
# Last modified: 5 Oct 2016

# Modified 5 Oct: added regional o-stats function.

# Takes as input an O-stats object.


ostat2longform <- function(o) {
  result_names <- c('site','trait','ostat_norm','ostat_norm_localnull_lower','ostat_norm_localnull_upper','ostat_norm_regnull_lower','ostat_norm_regnull_upper','ostat_norm_localnull_ses','ostat_norm_localnull_seslower','ostat_norm_localnull_sesupper','ostat_norm_regnull_ses','ostat_norm_regnull_seslower','ostat_norm_regnull_sesupper','ostat_unnorm','ostat_unnorm_localnull_lower','ostat_unnorm_localnull_upper','ostat_unnorm_regnull_lower','ostat_unnorm_regnull_upper','ostat_unnorm_localnull_ses','ostat_unnorm_localnull_seslower','ostat_unnorm_localnull_sesupper','ostat_unnorm_regnull_ses','ostat_unnorm_regnull_seslower','ostat_unnorm_regnull_sesupper')
  
  res_list <- list()
  
  nsite <- nrow(o[[1]])
  ntrait <- ncol(o[[1]])
  
  for (i in 1:nsite) {
    for (j in 1:ntrait) {
      res_list[[length(res_list)+1]] <- c(o$overlaps_norm[i,j],
                                          o$overlaps_norm_ses$raw_lower[i,j],
                                          o$overlaps_norm_ses$raw_upper[i,j],
                                          o$regional_overlaps_norm_ses$raw_lower[i,j],
                                          o$regional_overlaps_norm_ses$raw_upper[i,j],
                                          o$overlaps_norm_ses$ses[i,j],
                                          o$overlaps_norm_ses$ses_lower[i,j],
                                          o$overlaps_norm_ses$ses_upper[i,j],
                                          o$regional_overlaps_norm_ses$ses[i,j],
                                          o$regional_overlaps_norm_ses$ses_lower[i,j],
                                          o$regional_overlaps_norm_ses$ses_upper[i,j],
                                          o$overlaps_unnorm[i,j],
                                          o$overlaps_unnorm_ses$raw_lower[i,j],
                                          o$overlaps_unnorm_ses$raw_upper[i,j],
                                          o$regional_overlaps_unnorm_ses$raw_lower[i,j],
                                          o$regional_overlaps_unnorm_ses$raw_upper[i,j],
                                          o$overlaps_unnorm_ses$ses[i,j],
                                          o$overlaps_unnorm_ses$ses_lower[i,j],
                                          o$overlaps_unnorm_ses$ses_upper[i,j],
                                          o$regional_overlaps_unnorm_ses$ses[i,j],
                                          o$regional_overlaps_unnorm_ses$ses_lower[i,j],
                                          o$regional_overlaps_unnorm_ses$ses_upper[i,j])
    }
  }
  res <- as.data.frame(do.call('rbind', res_list))
  res <- cbind(site = rep(dimnames(o[[1]])[[1]], each=ntrait), trait = rep(dimnames(o[[1]])[[2]], times=nsite), res)
  names(res) <- result_names
  res
}


regionalostat2longform <- function(o) {
  result_names <- c('site','trait','ostat_reg','ostat_reg_allpoolnull_lower','ostat_reg_allpoolnull_upper','ostat_reg_byspnull_lower','ostat_reg_byspnull_upper','ostat_reg_allpoolnull_ses','ostat_reg_allpoolnull_seslower','ostat_reg_allpoolnull_sesupper','ostat_reg_byspnull_ses','ostat_reg_byspnull_seslower','ostat_reg_byspnull_sesupper')
  
  res_list <- list()
  
  nsite <- nrow(o[[1]])
  ntrait <- ncol(o[[1]])
  
  for (i in 1:nsite) {
    for (j in 1:ntrait) {
      res_list[[length(res_list)+1]] <- c(o$overlaps_reg[i,j],
                                          o$overlaps_reg_allpool_ses$raw_lower[i,j],
                                          o$overlaps_reg_allpool_ses$raw_upper[i,j],
                                          o$overlaps_reg_bysp_ses$raw_lower[i,j],
                                          o$overlaps_reg_bysp_ses$raw_upper[i,j],
                                          o$overlaps_reg_allpool_ses$ses[i,j],
                                          o$overlaps_reg_allpool_ses$ses_lower[i,j],
                                          o$overlaps_reg_allpool_ses$ses_upper[i,j],
                                          o$overlaps_reg_bysp_ses$ses[i,j],
                                          o$overlaps_reg_bysp_ses$ses_lower[i,j],
                                          o$overlaps_reg_bysp_ses$ses_upper[i,j]
                                          )
    }
  }
  res <- as.data.frame(do.call('rbind', res_list))
  res <- cbind(site = rep(dimnames(o[[1]])[[1]], each=ntrait), trait = rep(dimnames(o[[1]])[[2]], times=nsite), res)
  names(res) <- result_names
  res
}