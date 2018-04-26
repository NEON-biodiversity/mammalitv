# Variance partitioning functions
# Uses equations from Siefert et al. 2015, Ecology Letters, Box 1
# Written by QDR, March 2016


# Calculate ratio of intraspecific to total variance within communities.

wITV <- function(species, trait_values, abundance) {
  # Get rid of NA values and species with no data.
  abundance[is.na(abundance)] <- 0
  na_idx <- !is.na(trait_values)
  species <- as.character(species[na_idx])
  trait_values <- trait_values[na_idx]
  abundance <- abundance[na_idx]
  
  # Calculate intraspecific variance
  t_ij <- tapply(trait_values, species, mean)
  N_ij <- tapply(trait_values, species, function(x) sum(!is.na(x)))
  p_ij <- tapply(abundance, species, '[', 1)
  N_ij[is.na(N_ij)] <- 0 # Deals with species that have no recorded data.
  p_ij[is.na(p_ij)] <- 0
  ss_tijk <- sum((trait_values - t_ij[match(species, names(t_ij))])^2)
  intra_var <- sum(p_ij * (1/N_ij) * ss_tijk)
  
  # Calculate total variance
  total_var <- sum(p_ij * (t_ij - sum(p_ij*t_ij))^2) + intra_var
  return(c(wITV=intra_var/total_var, intra_var=intra_var, total_var=total_var))
}

# Calculate ratio of intraspecific to total among communities
aITV <- function(species, trait_values, abundance, sites) {
  # Get rid of NA values.
  abundance[is.na(abundance)] <- 0
  na_idx <- !is.na(trait_values)
  species <- species[na_idx]
  trait_values <- trait_values[na_idx]
  abundance <- abundance[na_idx]
  sites <- sites[na_idx]
  
  fixed_traitmean <- tapply(trait_values, species, mean)
  
  require(plyr)
  dat <- data.frame(species,trait_values,abundance,sites)
  dat_mean <- ddply(dat, .(sites, species, abundance), summarize, mean_i = mean(trait_values))
  dat_mean$mean_fixed <- fixed_traitmean[match(dat_mean$species, names(fixed_traitmean))]
  CWMs <- ddply(dat_mean, .(sites), summarize, CWM_i = weighted.mean(x=mean_i, w=abundance), CWM_fixed_i = weighted.mean(x=mean_fixed, w=abundance))
  CWM_i <- CWMs$CWM_i
  CWM_fixed_i <- CWMs$CWM_fixed_i
  
  CWM_intra_i <- CWM_i - CWM_fixed_i
  SS_tot <- anova(lm(CWM_i ~ 1))$Sum
  SS_intra <- anova(lm(CWM_intra_i ~ 1))$Sum
  SS_fixed <- anova(lm(CWM_fixed_i ~ 1))$Sum
  return(c(aITV = log(SS_intra/SS_fixed), cov = 100 * (SS_tot-SS_intra-SS_fixed)/SS_tot, SS_tot=SS_tot, SS_intra=SS_intra, SS_fixed=SS_fixed))
}