# Code for archiving with Proc B manuscript
# 24 July 2017 QDR

# Tested under R/3.3.3

# Necessary steps to replicate analysis in MS
# -------------------------------------------

# 1. Define functions -----------------------------------------------------

# Overlap between two empirical density estimates
# Link: http://stats.stackexchange.com/questions/97596/how-to-calculate-overlap-between-empirical-probability-densities
pairwise_overlap <- function(a, b, norm = TRUE, bw = NULL, n = NULL) {
  
  # clean input
  a <- as.numeric(na.omit(a))
  b <- as.numeric(na.omit(b))
  
  # define limits of a common grid, adding a buffer so that tails aren't cut off
  lower <- min(c(a, b)) - 1 
  upper <- max(c(a, b)) + 1
  
  # generate kernel densities
  # add option to use user-defined bandwidth and n
  if (is.null(bw)) bw <- 'nrd0' # Defaults to traditional method if not given
  if (is.null(n)) n <- 512 # Default value if not given
  da <- density(a, from=lower, to=upper, bw=bw, n=n)
  db <- density(b, from=lower, to=upper, bw=bw, n=n)
  d <- data.frame(x=da$x, a=da$y, b=db$y)
  
  # If not normalized, multiply each density entry by the length of each vector
  if (!norm) {
    d$a <- d$a * length(a)
    d$b <- d$b * length(b)
  }
  
  # calculate intersection densities
  d$w <- pmin(d$a, d$b)
  
  # integrate areas under curves
  suppressMessages(require(sfsmisc))
  total <- integrate.xy(d$x, d$a) + integrate.xy(d$x, d$b)
  intersection <- integrate.xy(d$x, d$w)
  
  # compute overlap coefficient
  overlap <- 2 * intersection / total
  overlap_a <- intersection / integrate.xy(d$x, d$a)
  overlap_b <- intersection / integrate.xy(d$x, d$b)
  
  return(c(overlap = overlap, overlap_a = overlap_a, overlap_b = overlap_b))
  
}

# Median pairwise overlap of trait distributions of all species in a community
# The median is weighted by the harmonic mean of the abundances of each species pair for which overlap is calculated.
community_overlap_harmonicwmedian <- function(traits, sp, norm = TRUE, bw = NULL, n = NULL, randomize_weights = FALSE) {
  require(matrixStats)
  sp <- as.character(sp)
  dat <- data.frame(traits=traits, sp=sp, stringsAsFactors = FALSE)
  dat <- dat[complete.cases(dat), ]
  abunds <- table(dat$sp)
  abunds <- abunds[abunds>1]
  dat <- dat[dat$sp %in% names(abunds), ]
  traitlist <- split(dat$traits, dat$sp)
  nspp <- length(traitlist)
  
  if (nspp < 2) return(NA)
  
  overlaps <- numeric(0)
  abund_pairs <- numeric(0)
  
  for (sp_a in 1:(nspp-1)) {
    for (sp_b in (sp_a+1):nspp) {
      o <- pairwise_overlap(a = traitlist[[sp_a]], b = traitlist[[sp_b]], norm = norm, bw = bw, n = n)
      overlaps <- c(overlaps, o[1])
      harmonic_mean <- 2/(1/abunds[sp_a] + 1/abunds[sp_b])
      abund_pairs <- c(abund_pairs, harmonic_mean)
    }
  }
  
  if (randomize_weights) abund_pairs <- sample(abund_pairs)
  
  weightedMedian(x = overlaps, w = abund_pairs)
  
}

# Calculation of standardized effect sizes from null model values
get_ses <- function(o, o_n, qs) {
  ses <- ses_lower <- ses_upper <- raw_lower <- raw_upper <- matrix(NA, nrow=nrow(o), ncol=ncol(o))
  
  for (i in 1:nrow(o)) {
    for (j in 1:ncol(o)) {
      if(!is.na(o[i,j])) {
        obs <- o[i,j]
        nullvals <- na.omit(o_n[i, j, ])
        ses[i,j] <- (obs - mean(nullvals))/sd(nullvals)
        ses_lower[i,j] <- (quantile(nullvals, probs=qs[1]) - mean(nullvals))/sd(nullvals)
        ses_upper[i,j] <- (quantile(nullvals, probs=qs[2]) - mean(nullvals))/sd(nullvals)
        raw_lower[i,j] <- quantile(nullvals, probs=qs[1])
        raw_upper[i,j] <- quantile(nullvals, probs=qs[2])
      }
    }
  }
  dimnames(ses) <- dimnames(ses_lower) <- dimnames(ses_upper) <- dimnames(raw_lower) <- dimnames(raw_upper) <- dimnames(o)
  return(list(ses=ses, ses_lower=ses_lower, ses_upper=ses_upper, raw_lower=raw_lower, raw_upper=raw_upper))
}

# Calculation of overlap statistics, and implementation of null model that swaps means of distributions 
# to test whether distributions are more evenly spaced than expected by chance

Ostats <- function(traits, plots, sp, nperm = 99, nullqs = c(0.025, 0.975), shuffle_weights = FALSE, swap_means = FALSE) {
  # Required input: a matrix called traits (nrows=n individuals, ncols=n traits), 
  # a vector called plots which is a factor with length equal to nrow(traits),
  # a vector called sp which is a factor with length equal to nrow(traits),

  # Declaration of data structures to hold the results
  
  # Data structures for observed O-Stats
  overlaps_norm <- matrix(nrow = nlevels(plots), ncol = ncol(traits))
  overlaps_unnorm <- matrix(nrow = nlevels(plots), ncol = ncol(traits))

  # Data structures for null O-Stats
  overlaps_norm_null <- array(dim = c(nlevels(plots), ncol(traits), nperm))
  overlaps_unnorm_null <- array(dim = c(nlevels(plots), ncol(traits), nperm))

  # Name rows and columns of the outputs.
  dimnames(overlaps_norm) <- list(levels(plots), dimnames(traits)[[2]])
  dimnames(overlaps_unnorm) <- list(levels(plots), dimnames(traits)[[2]])

  # Calculation of observed O-Stats
  
  print('Calculating observed local O-stats for each community . . .')
  pb <- txtProgressBar(min = 0, max = nlevels(plots), style = 3)
  
  for (s in 1:nlevels(plots)) {
    for (t in 1:ncol(traits)) {
      overlap_norm_st <- try(community_overlap_harmonicwmedian(traits = traits[plots == levels(plots)[s], t], sp = sp[plots == levels(plots)[s]], norm=TRUE), TRUE)
      overlaps_norm[s, t] <- if (inherits(overlap_norm_st, 'try-error')) NA else overlap_norm_st
      overlap_unnorm_st <- try(community_overlap_harmonicwmedian(traits = traits[plots == levels(plots)[s], t], sp = sp[plots == levels(plots)[s]], norm=FALSE), TRUE)
      overlaps_unnorm[s, t] <- if (inherits(overlap_unnorm_st, 'try-error')) NA else overlap_unnorm_st
    }
    setTxtProgressBar(pb, s)
  }
  
  close(pb)

  print('Calculating null distributions of O-stats . . . ')
  pb <- txtProgressBar(min = 0, max = nperm, style = 3)
  
  # Null model generation and calculation of null O-Stats
  
  # Local null model: generation and calculation done in the same loop
  for (i in 1:nperm) {
    setTxtProgressBar(pb, i)
    for (s in 1:nlevels(plots)) {
      for (t in 1:ncol(traits)) {
          if (!shuffle_weights & !swap_means) overlap_norm_sti <- try(community_overlap_harmonicwmedian(traits = traits[plots == levels(plots)[s], t], sp = sample(sp[plots == levels(plots)[s]]), norm=TRUE), TRUE)
          if (shuffle_weights) overlap_norm_sti <- try(community_overlap_harmonicwmedian(traits = traits[plots == levels(plots)[s], t], sp = sp[plots == levels(plots)[s]], norm=TRUE, randomize_weights = TRUE), TRUE)
          if (swap_means) {
            traits_st <- traits[plots==levels(plots)[s], t]
            sp_st <- sp[plots==levels(plots)[s]]
            
            traitmeans <- tapply(traits_st,sp_st,mean)
            traitdeviations <- traits_st-traitmeans[sp_st]
            
            # Sort the trait means out randomly.
            traitmeans_null <- sample(traitmeans)
            sp_null <- rep(names(traitmeans_null), table(sp_st))
            traits_null <- traitdeviations + traitmeans_null[sp_null]
            overlap_norm_sti <- try(community_overlap_harmonicwmedian(traits = traits_null, sp = sp_null, norm=TRUE, randomize_weights = FALSE), TRUE)
          }
        
        overlaps_norm_null[s, t, i] <- if (inherits(overlap_norm_sti, 'try-error')) NA else overlap_norm_sti
        overlap_unnorm_sti <- try(community_overlap_harmonicwmedian(traits = traits[plots == levels(plots)[s], t], sp = sample(sp[plots == levels(plots)[s]]), norm=FALSE), TRUE)
        overlaps_unnorm_null[s, t, i] <- if (inherits(overlap_unnorm_sti, 'try-error')) NA else overlap_unnorm_sti
      }
    }
  }  
  
  close(pb) 
  print('Extracting null quantiles to get standardized effect sizes (almost done!) . . .')
  
  # Extract quantiles to get standardized effect sizes for the overlap stats
  overlaps_norm_ses <- get_ses(overlaps_norm, overlaps_norm_null, nullqs)
  overlaps_unnorm_ses <- get_ses(overlaps_unnorm, overlaps_unnorm_null, nullqs)
  list(overlaps_norm=overlaps_norm, overlaps_unnorm=overlaps_unnorm, 
       overlaps_norm_ses=overlaps_norm_ses, overlaps_unnorm_ses=overlaps_unnorm_ses)
  
}

# Convert output of overlap statistics to long form
ostat2longform <- function(o) {
  result_names <- c('site','trait','ostat_norm','ostat_norm_localnull_lower','ostat_norm_localnull_upper','ostat_norm_localnull_ses','ostat_norm_localnull_seslower','ostat_norm_localnull_sesupper','ostat_unnorm','ostat_unnorm_localnull_lower','ostat_unnorm_localnull_upper','ostat_unnorm_localnull_ses','ostat_unnorm_localnull_seslower','ostat_unnorm_localnull_sesupper')
  
  res_list <- list()
  
  nsite <- nrow(o[[1]])
  ntrait <- ncol(o[[1]])
  
  for (i in 1:nsite) {
    for (j in 1:ntrait) {
      res_list[[length(res_list)+1]] <- c(o$overlaps_norm[i,j],
                                          o$overlaps_norm_ses$raw_lower[i,j],
                                          o$overlaps_norm_ses$raw_upper[i,j],
                                          o$overlaps_norm_ses$ses[i,j],
                                          o$overlaps_norm_ses$ses_lower[i,j],
                                          o$overlaps_norm_ses$ses_upper[i,j],
                                          o$overlaps_unnorm[i,j],
                                          o$overlaps_unnorm_ses$raw_lower[i,j],
                                          o$overlaps_unnorm_ses$raw_upper[i,j],
                                          o$overlaps_unnorm_ses$ses[i,j],
                                          o$overlaps_unnorm_ses$ses_lower[i,j],
                                          o$overlaps_unnorm_ses$ses_upper[i,j])
    }
  }
  res <- as.data.frame(do.call('rbind', res_list))
  res <- cbind(site = rep(dimnames(o[[1]])[[1]], each=ntrait), trait = rep(dimnames(o[[1]])[[2]], times=nsite), res)
  names(res) <- result_names
  res
}


# 2. Load necessary packages ----------------------------------------------

library(dplyr) # Data manipulation
library(lubridate) # Data manipulation
library(iNEXT) # Calculation of richness estimators
library(MuMIn) # Model selection
library(matrixStats) # Calculation of weighted median
library(blavaan) # Fitting structural equation models

# 3. Load data ------------------------------------------------------------

# Climate and landscape covariates at site level, averaged from all plot-level covariates, including pre-calculated PCA axes
site_covariates <- read.csv('site_covariates.csv', stringsAsFactors = FALSE) 

# Mammal trapping data from NEON, including all 2015 data
raw_mammal_data <- read.csv('raw_NEON_mammal_data.csv', stringsAsFactors = FALSE)

# Quality-controlled NEON data
final_mammal_data <- read.csv('final_NEON_mammal_data.csv', stringsAsFactors = FALSE)

# Mammal trait data compiled from various external sources
mammal_traits <- read.csv('mammal_traits.csv', stringsAsFactors = FALSE)

# 4. Quality control on NEON data -----------------------------------------

sites_to_exclude <- c('DSNY', 'DELA', 'HEAL') # See manuscript for why these sites are excluded

# Get rid of all non-target taxa and combine different columns that were used for individual IDs in different contexts.
mammal_data <- mutate(raw_mammal_data, 
                      individualandtag = pmin(as.character(individualID), as.character(tagID), na.rm=TRUE),
                      year = year(date)) %>%
  filter(order == 'Rodentia', taxonProtocolCategory == 'target', year == 2015)


# 5. Calculate richness estimators ----------------------------------------

# Get rid of poorly identified individuals (those marked sp.) because they should not count for our sampling . . . they are spurious singletons.
# Then generate vectors of abundances by species for each site.
mammaltables <- mammal_data %>% 
  filter(!siteID %in% sites_to_exclude, !grepl('sp\\.', mammal_data$scientificName)) %>% 
  group_by(siteID) %>% 
  do(t = table(.$taxonID))

# Name the list of vectors
mamx <- lapply(mammaltables$t, as.numeric)
names(mamx) <- mammaltables$siteID

# Calculate asymptotic richness estimator
set.seed(46545)
richness_estimators <- iNEXT(x=mamx, q=0, datatype='abundance', size = c(5,10,50,100,500,1000,2000,3000,4000))

# Calculate Chao1 richness estimator and combine all richness estimators by site.
chao <- function(x) {
  xcomm <- table(x$taxonID)
  S_obs <- length(xcomm)
  f1 <- sum(xcomm == 1)
  f2 <- sum(xcomm == 2)
  return(data.frame(chao1 = S_obs + (f1 * (f1 - 1)) / (2 * (f2 + 1))))
}

asymptotic_richness <- richness_estimators$AsyEst %>% filter(Diversity == 'Species richness') 
asymptotic_richness$Site <- factor(asymptotic_richness$Site, levels=asymptotic_richness$Site[order(asymptotic_richness$Observed)])

chao1site <- mammal_data %>% 
  filter(!siteID %in% sites_to_exclude, !grepl('sp\\.', scientificName)) %>% 
  group_by(siteID) %>% 
  do(chao(.))

chao1site$siteID <-factor(asymptotic_richness$Site, levels=asymptotic_richness$Site[order(asymptotic_richness$Observed)])

richnessdat <- full_join(asymptotic_richness %>% select(Site, Observed, Estimator) %>% rename(siteID=Site), chao1site)



# 6. Calculate overlap statistics and null effect sizes -------------------

mammal_logweight <- final_mammal_data[,c('logweight'), drop = FALSE]

set.seed(37917)
Ostats_bysite2015 <- with(final_mammal_data,
                          Ostats(traits = mammal_logweight, plots = factor(siteID), sp = factor(taxonID),  nperm = 999, shuffle_weights = FALSE, swap_means = TRUE))

ostats <- ostat2longform(Ostats_bysite2015)
ostats <- ostats %>% 
  filter(!site %in% sites_to_exclude) %>% 
  rename(siteID=site) %>% 
  left_join(site_covariates) %>% 
  left_join(richnessdat) 


# 7. Fit SEMs -------------------------------------------------------------

fullsemdat <- ostats %>%
  mutate(o = qlogis(ostat_norm),
         r = (chao1 - mean(chao1, na.rm=T))/sd(chao1, na.rm=T),
         p = (log10(NPP) - mean(log10(NPP), na.rm=T))/sd(log10(NPP), na.rm=T),
         tmin = (bio6 - mean(bio6, na.rm=T))/sd(bio6, na.rm=T)) %>%
 select(siteID, o, r, p, tmin) 

# Model specifications for SEMs

# This is the full model with indirect (overlap-mediated) effects and direct effects for both temperature (bio6) and productivity (NPP)
tempdi_proddi <- 
' # direct effects
r ~ e*tmin
r ~ f*p
# mediator
o ~ a*tmin + b*p
r ~ c*o
# indirect effect of temperature (a*c)
indir_temp := a*c
# total effect of temperature
total_temp := e + (a*c)
# indirect effect of productivity (b*c)
indir_prod := b*c
# total effect of productivity
total_prod := f + (b*c)
# intercepts
o ~ 1
r ~ 1
'

# Indirect and direct temperature effects, with no productivity effects
tempdi <- 
' # direct effects
r ~ e*tmin
# mediator
o ~ a*tmin
r ~ c*o
# indirect effect of temperature (a*c)
indir_temp := a*c
# total effect of temperature
total_temp := e + (a*c)
# intercepts
o ~ 1
r ~ 1
'

# Indirect temperature effect, with no productivity effects
tempi <- 
' 
# mediator
o ~ a*tmin
r ~ c*o
# indirect effect of temperature (a*c)
indir_temp := a*c
# intercepts
o ~ 1
r ~ 1
'

# Direct temperature effect, with no productivity effects
tempd <- 
' # direct effects
r ~ e*tmin
# mediator
r ~ c*o
# intercepts
r ~ 1
'

# Direct and indirect temperature effects, only indirect productivity effect
tempdi_prodi <- 
' # direct effects
r ~ e*tmin
# mediator
o ~ a*tmin + b*p
r ~ c*o
# indirect effect of temperature (a*c)
indir_temp := a*c
# total effect of temperature
total_temp := e + (a*c)
# indirect effect of productivity (b*c)
indir_prod := b*c
# intercepts
o ~ 1
r ~ 1
'

# Indirect temperature effect, direct and indirect productivity effects
tempi_proddi <- 
' # direct effects
r ~ f*p
# mediator
o ~ a*tmin + b*p
r ~ c*o
# indirect effect of temperature (a*c)
indir_temp := a*c
# indirect effect of productivity (b*c)
indir_prod := b*c
# total effect of productivity
total_prod := f + (b*c)
# intercepts
o ~ 1
r ~ 1
'

# Indirect temperature effect, indirect productivity effect
tempi_prodi <- 
' 
# mediator
o ~ a*tmin + b*p
r ~ c*o
# indirect effect of temperature (a*c)
indir_temp := a*c
# indirect effect of productivity (b*c)
indir_prod := b*c
# intercepts
o ~ 1
r ~ 1
'

# Run SEM model fits

set.seed(4003846)

tempdi_proddi_fit <- blavaan(model = tempdi_proddi, data = fullsemdat,
                             auto.var = TRUE, auto.fix.first = TRUE, auto.cov.lv.x = TRUE,
                             jagcontrol = list(method = 'rjparallel'),
                             n.chains = 3, burnin = 5000, sample = 20000)
tempdi_fit <- blavaan(model = tempdi, data = fullsemdat,
                      auto.var = TRUE, auto.fix.first = TRUE, auto.cov.lv.x = TRUE,
                      jagcontrol = list(method = 'rjparallel'),
                      n.chains = 3, burnin = 5000, sample = 20000)
tempi_fit <- blavaan(model = tempi, data = fullsemdat,
                     auto.var = TRUE, auto.fix.first = TRUE, auto.cov.lv.x = TRUE,
                     jagcontrol = list(method = 'rjparallel'),
                     n.chains = 3, burnin = 5000, sample = 20000)
tempd_fit <- blavaan(model = tempd, data = fullsemdat,
                     auto.var = TRUE, auto.fix.first = TRUE, auto.cov.lv.x = TRUE,
                     jagcontrol = list(method = 'rjparallel'),
                     n.chains = 3, burnin = 5000, sample = 20000)
tempdi_prodi_fit <- blavaan(model = tempdi_prodi, data = fullsemdat,
                            auto.var = TRUE, auto.fix.first = TRUE, auto.cov.lv.x = TRUE,
                            jagcontrol = list(method = 'rjparallel'),
                            n.chains = 3, burnin = 5000, sample = 20000)
tempi_proddi_fit <- blavaan(model = tempi_proddi, data = fullsemdat,
                            auto.var = TRUE, auto.fix.first = TRUE, auto.cov.lv.x = TRUE,
                            jagcontrol = list(method = 'rjparallel'),
                            n.chains = 3, burnin = 5000, sample = 20000)
tempi_prodi_fit <- blavaan(model = tempi_prodi, data = fullsemdat,
                           auto.var = TRUE, auto.fix.first = TRUE, auto.cov.lv.x = TRUE,
                           jagcontrol = list(method = 'rjparallel'),
                           n.chains = 3, burnin = 5000, sample = 20000)

# 8. Extract fit statistics and compare SEMs ------------------------------

all_sem_stats <- function(fit) {
  cat('Parameter estimates\n-------------------\n')
  print(parameterEstimates(fit))
  cat('\nR-square\n--------\n')
  print(inspect(fit, 'rsquare'))
  cat('\nFit measures\n------------\n')
  print(fitMeasures(fit))
  cat('\nPseudo-p-values\n---------------\n')
  
  allchains <- do.call('rbind', fit@external$runjags$mcmc)
  print(apply(allchains>0,2,sum)/nrow(allchains))
  print(apply(allchains<0,2,sum)/nrow(allchains))
  
}

all_sem_stats(tempdi_proddi_fit)
all_sem_stats(tempdi_fit)
all_sem_stats(tempi_fit)
all_sem_stats(tempd_fit)
all_sem_stats(tempdi_prodi_fit)
all_sem_stats(tempi_proddi_fit)
all_sem_stats(tempi_prodi_fit)

# 9. Calculate community- and species-level quantile ranges ---------------

qprobs <- c(0.05, 0.95)

commiqrs <- final_mammal_data %>% 
  transform(logweight = log10(weight)) %>%
  group_by(siteID) %>%
  summarize(community_iqr = diff(quantile(logweight, probs = qprobs, na.rm = T)))

spiqrs <- final_mammal_data %>% 
  transform(logweight = log10(weight)) %>%
  group_by(siteID, taxonID) %>%
  summarize(species_iqr = diff(quantile(logweight, probs = qprobs, na.rm = T)),
            species_n = sum(!is.na(logweight)))

spiqrs_mean <- spiqrs %>% ungroup %>%
  group_by(siteID) %>%
  summarize(byspecies_iqr = weighted.mean(species_iqr, species_n, na.rm=T))

iqrdat <- ostats %>%
  dplyr::select(siteID, chao1, bio1, bio6, ostat_norm) %>%
  left_join(commiqrs) %>% left_join(spiqrs_mean)

iqrdat$niche_zscore <- ostats$ostat_norm_localnull_ses

# 10. Run linear model for three drivers of overlap -----------------------

lm_maineffects_standard <- lm(I(qlogis(ostat_norm)) ~ scale(community_iqr) + scale(byspecies_iqr) + scale(niche_zscore), data=iqrdat %>% filter(complete.cases(.)), na.action='na.pass')

dredge(lm_maineffects_standard)
summary(lm_maineffects_standard)
confint(lm_maineffects_standard)

