################################################################################################################
################################################################################################################
################################################################################################################
#### R code to reproduce all analyses in the paper:							                            	####
#### Intraspecific variation in body size reveals community assembly in North American rodents              ####
################################################################################################################
#### Quentin D. Read, John M. Grady, Phoebe Zarnetske, Sydne Record, Benjamin Baiser, Jonathan Belmaker,    ####
#### Mao-Ning Tuanmu, Angela Strecker, Lydia Beaudrot, Kate Thibault                       					####
################################################################################################################
#### Contact: Quentin Read, qdr@msu.edu									              						####
################################################################################################################
#### Last modified: 16 January 2017           		                                    					####
################################################################################################################
################################################################################################################
################################################################################################################

# This code was tested with R version 3.3.1 on 16 January 2017.
# Required packages: dplyr, sfsmisc, matrixStats, betareg, MuMIn, usdm, iNEXT, lavaan

# Contents of this document:
#	I. Functions needed to analyze data
#	II. Input and processing of community trait data and site-level covariate data
#	III. Calculation of overlap statistics and generation of null model quantiles.
#	IV. Regression of overlap statistics on environmental variables.
#	V. Fit alternative path models relating temperature, richness, and overlap.

##############################################	
## Part I: Functions needed to analyze data ##
##############################################

# 1. Overlap of two empirical probability densities
# Adapted from code on stack overflow
# Link: http://stats.stackexchange.com/questions/97596/how-to-calculate-overlap-between-empirical-probability-densities

pairwise_overlap <- function(a, b, norm=TRUE, bw = NULL, n = NULL) {
  
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

# 2. Calculate abundance-weighted community-level median of pairwise trait overlap (abundance weights are the harmonic means of the abundances of each species pair)

community_overlap_harmonicwmedian <- function(traits, sp, norm = TRUE, bw = NULL, n = NULL, randomize_weights = FALSE) {
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
  
  matrixStats::weightedMedian(x = overlaps, w = abund_pairs)
  
}


# 3. Function to extract null model quantiles from the generated distributions of overlap statistics.

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

# 4. Function to calculate overlap statistics for a series of communities
# Note: there are a number of additional statistics output by this function that do not appear in the manuscript.
# The metrics based on the regional species pool may be used in future analyses of the NEON data.

Ostats <- function(traits, plots, sp, reg_pool_traits, reg_pool_sp, nperm = 99, nullqs = c(0.025, 0.975), stat = 'median', shuffle_weights = FALSE) {
	# Required input: a matrix called traits (nrows=n individuals, ncols=n traits), 
	# a vector called plots which is a factor with length equal to nrow(traits),
	# a vector called sp which is a factor with length equal to nrow(traits),
	# a list called reg_pool_traits of which each element is a matrix with ncols=n traits
	# a list called reg_pool_sp of which each element is a factor with length equal to the 
	# number of rows in the corresponding reg_pool_traits matrix.

	# Declaration of data structures to hold the results
	
	# Data structures for observed O-Stats
	overlaps_norm <- matrix(nrow = nlevels(plots), ncol = ncol(traits))
	overlaps_unnorm <- matrix(nrow = nlevels(plots), ncol = ncol(traits))
	#regional_overlaps_norm <- matrix(nrow = nlevels(plots), ncol = ncol(traits))
	#regional_overlaps_unnorm <- matrix(nrow = nlevels(plots), ncol = ncol(traits))
	
	# Data structures for null O-Stats
	overlaps_norm_null <- array(dim = c(nlevels(plots), ncol(traits), nperm))
	overlaps_unnorm_null <- array(dim = c(nlevels(plots), ncol(traits), nperm))
	regional_overlaps_norm_null <- array(dim = c(nlevels(plots), ncol(traits), nperm))
	regional_overlaps_unnorm_null <- array(dim = c(nlevels(plots), ncol(traits), nperm))
	
	# Name rows and columns of the outputs.
	dimnames(overlaps_norm) <- list(levels(plots), dimnames(traits)[[2]])
	dimnames(overlaps_unnorm) <- list(levels(plots), dimnames(traits)[[2]])
	
	# Calculation of observed O-Stats
	
	print('Calculating observed local O-stats for each community . . .')
	pb <- txtProgressBar(min = 0, max = nlevels(plots), style = 3)
	
	for (s in 1:nlevels(plots)) {
		for (t in 1:ncol(traits)) {
			if (stat == 'median') overlap_norm_st <- try(community_overlap(traits = traits[plots == levels(plots)[s], t], sp = sp[plots == levels(plots)[s]], norm=TRUE), TRUE)
			if (stat == 'weighted.mean') overlap_norm_st <- try(community_overlap_wm(traits = traits[plots == levels(plots)[s], t], sp = sp[plots == levels(plots)[s]], norm=TRUE), TRUE)
			if (stat == 'weighted.median') overlap_norm_st <- try(community_overlap_wmedian(traits = traits[plots == levels(plots)[s], t], sp = sp[plots == levels(plots)[s]], norm=TRUE), TRUE)
			if (stat == 'harmonic') overlap_norm_st <- try(community_overlap_harmonicwmedian(traits = traits[plots == levels(plots)[s], t], sp = sp[plots == levels(plots)[s]], norm=TRUE), TRUE)
			if (stat == 'none') overlap_norm_st <- try(community_overlap_noitv(traits = traits[plots == levels(plots)[s], t], sp = sp[plots == levels(plots)[s]], norm=TRUE), TRUE)
			overlaps_norm[s, t] <- if (inherits(overlap_norm_st, 'try-error')) NA else overlap_norm_st
			overlap_unnorm_st <- try(community_overlap(traits = traits[plots == levels(plots)[s], t], sp = sp[plots == levels(plots)[s]], norm=FALSE), TRUE)
			overlaps_unnorm[s, t] <- if (inherits(overlap_unnorm_st, 'try-error')) NA else overlap_unnorm_st
		}
		setTxtProgressBar(pb, s)
	}
	
	close(pb)

	print('Calculating local and regional null distributions of O-stats . . . ')
	pb <- txtProgressBar(min = 0, max = nperm, style = 3)
	
	# Null model generation and calculation of null O-Stats
	
	# Local null model: generation and calculation done in the same loop
	for (i in 1:nperm) {
		for (s in 1:nlevels(plots)) {
			for (t in 1:ncol(traits)) {
				if (stat == 'median') overlap_norm_sti <- try(community_overlap(traits = traits[plots == levels(plots)[s], t], sp = sample(sp[plots == levels(plots)[s]]), norm=TRUE), TRUE)
				if (stat == 'weighted.mean') overlap_norm_sti <- try(community_overlap_wm(traits = traits[plots == levels(plots)[s], t], sp = sample(sp[plots == levels(plots)[s]]), norm=TRUE), TRUE)
				if (stat == 'weighted.median') overlap_norm_sti <- try(community_overlap_wmedian(traits = traits[plots == levels(plots)[s], t], sp = sample(sp[plots == levels(plots)[s]]), norm=TRUE), TRUE)
				if (stat == 'harmonic') {
					if (!shuffle_weights) overlap_norm_sti <- try(community_overlap_harmonicwmedian(traits = traits[plots == levels(plots)[s], t], sp = sample(sp[plots == levels(plots)[s]]), norm=TRUE), TRUE)
					if (shuffle_weights) overlap_norm_sti <- try(community_overlap_harmonicwmedian(traits = traits[plots == levels(plots)[s], t], sp = sp[plots == levels(plots)[s]], norm=TRUE, randomize_weights = TRUE), TRUE)
				}
				if (stat == 'none') overlap_norm_sti <- try(community_overlap_noitv(traits = traits[plots == levels(plots)[s], t], sp = sample(sp[plots == levels(plots)[s]]), norm=TRUE), TRUE)
				overlaps_norm_null[s, t, i] <- if (inherits(overlap_norm_sti, 'try-error')) NA else overlap_norm_sti
				overlap_unnorm_sti <- try(community_overlap(traits = traits[plots == levels(plots)[s], t], sp = sample(sp[plots == levels(plots)[s]]), norm=FALSE), TRUE)
				overlaps_unnorm_null[s, t, i] <- if (inherits(overlap_unnorm_sti, 'try-error')) NA else overlap_unnorm_sti
			}
		}
		
	# Regional null model: generation of null communities and calculation of O-stats (in same loop)
	
	# Sample a number of individuals from the regional pool equal to the number of individuals in the community.		
				
		for (s in 1:nlevels(plots)) {
			for (t in 1:ncol(traits)) {
			
				n_indiv_st <- sum(plots == levels(plots)[s])		
				null_comm_index <- sample(x = length(reg_pool_sp[[s]]), size = n_indiv_st, replace = FALSE)
				
				if (stat == 'median') regional_overlap_norm_sti <- try(community_overlap(traits = reg_pool_traits[[s]][null_comm_index, t], sp = reg_pool_sp[[s]][null_comm_index], norm=TRUE), TRUE)
				if (stat == 'weighted.mean') regional_overlap_norm_sti <- try(community_overlap_wm(traits = reg_pool_traits[[s]][null_comm_index, t], sp = reg_pool_sp[[s]][null_comm_index], norm=TRUE), TRUE)
				if (stat == 'weighted.median') regional_overlap_norm_sti <- try(community_overlap_wmedian(traits = reg_pool_traits[[s]][null_comm_index, t], sp = reg_pool_sp[[s]][null_comm_index], norm=TRUE), TRUE)
				if (stat == 'harmonic') regional_overlap_norm_sti <- try(community_overlap_harmonicwmedian(traits = reg_pool_traits[[s]][null_comm_index, t], sp = reg_pool_sp[[s]][null_comm_index], norm=TRUE), TRUE)
				if (stat == 'none') regional_overlap_norm_sti <- try(community_overlap_noitv(traits = reg_pool_traits[[s]][null_comm_index, t], sp = reg_pool_sp[[s]][null_comm_index], norm=TRUE), TRUE)
				regional_overlaps_norm_null[s, t, i] <- if (inherits(regional_overlap_norm_sti, 'try-error')) NA else regional_overlap_norm_sti
				regional_overlap_unnorm_sti <- try(community_overlap(traits = reg_pool_traits[[s]][null_comm_index, t], sp = reg_pool_sp[[s]][null_comm_index], norm=FALSE), TRUE)
				regional_overlaps_unnorm_null[s, t, i] <- if (inherits(regional_overlap_unnorm_sti, 'try-error')) NA else regional_overlap_unnorm_sti
			}
		}
		setTxtProgressBar(pb, i)
	}
	
	close(pb)
	print('Extracting null quantiles to get standardized effect sizes (almost done!) . . .')
	
	# Extract quantiles to get standardized effect sizes for the overlap stats
		
	overlaps_norm_ses <- get_ses(overlaps_norm, overlaps_norm_null, nullqs)
	overlaps_unnorm_ses <- get_ses(overlaps_unnorm, overlaps_unnorm_null, nullqs)
	regional_overlaps_norm_ses <- get_ses(overlaps_norm, regional_overlaps_norm_null, nullqs)
	regional_overlaps_unnorm_ses <- get_ses(overlaps_unnorm, regional_overlaps_unnorm_null, nullqs)
	
	list(overlaps_norm=overlaps_norm, overlaps_unnorm=overlaps_unnorm, 
		 overlaps_norm_ses=overlaps_norm_ses, overlaps_unnorm_ses=overlaps_unnorm_ses, regional_overlaps_norm_ses=regional_overlaps_norm_ses, regional_overlaps_unnorm_ses=regional_overlaps_unnorm_ses)
	
}

# 5. Function to convert output from Ostats() to a long-form data frame

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

##########################################################################################
## Part II: Input and processing of community trait data and site-level covariate data. ##
##########################################################################################

library(dplyr)

# 1. Load site-level climate data
# This already contains z-scores of phylogenetic diversity, as well as the pre-calculated PCA axes for environmental heterogeneity.

siteclimatemeans <- read.csv('siteclimatemeans.csv', stringsAsFactors = FALSE)

# 2. Load community data

mammaltraits <- read.csv('mammaltraits.csv', stringsAsFactors = FALSE)

# Load raw data from NEON
mam_capture <- read.csv('mam_capture.csv', stringsAsFactors = FALSE)

# 3. Calculate richness estimators

# Chao1 richness estimator
chao <- function(x) {
  xcomm <- table(x$taxonID)
  S_obs <- length(xcomm)
  f1 <- sum(xcomm == 1)
  f2 <- sum(xcomm == 2)
  return(data.frame(chao1 = S_obs + (f1 * (f1 - 1)) / (2 * (f2 + 1))))
}

# Calculate both the asymptotic richness estimator and the Chao1 estimator.
# Exclude all non-target individuals and those that are only identified to the genus level.
library(iNEXT)

mammaltables <- mam_capture %>% filter(!siteID %in% c('DSNY','DELA','HEAL'), taxonProtocolCategory == 'target', order == 'Rodentia', !grepl('sp\\.', scientificName)) %>% group_by(siteID) %>% do(t = table(.$taxonID))

mamx <- lapply(mammaltables$t, as.numeric)
names(mamx) <- mammaltables$siteID

set.seed(46545)
default <- iNEXT(x=mamx, q=0, datatype='abundance', size = c(5,10,50,100,500,1000,2000,3000,4000))

asyrich <- default$AsyEst %>% filter(Diversity == 'Species richness') 
asyrich$Site <- factor(asyrich$Site, levels=asyrich$Site[order(asyrich$Observed)])

chao1site <- mam_capture %>% filter(!siteID %in% c('DSNY','DELA','HEAL'), taxonProtocolCategory == 'target', order == 'Rodentia', !grepl('sp\\.', scientificName)) %>% group_by(siteID) %>% do(chao(.))
chao1site$siteID <-factor(asyrich$Site, levels=asyrich$Site[order(asyrich$Observed)])

richness <- asyrich %>% select(Site, Observed, Estimator) %>% rename(siteID = Site, observed_richness = Observed, richness_estimator = Estimator) %>% left_join(chao1site)

#########################################################################################
## Part III: Calculation of overlap statistics and generation of null model quantiles. ##
#########################################################################################

# Set number of null model iterations here and which stat is used.
STAT <- 'harmonic'
N_PERM <- 999

# Process data frame for calculating statistics.
mammaltraits <- mammaltraits %>% filter(Pineda_Main_food != 'Insectivore' & !(siteID %in% c('DSNY','DELA','HEAL')))

# Retain only 2015 records
i2015 <- mammaltraits$year == 2015

# Convert mammal data to trait-only format
mammaltraits$logweight <- log10(mammaltraits$weight)
traits_mam15 <- mammaltraits[i2015, c('weight', 'logweight')] # Run on both raw and log-transformed weights.

sites15 <- unique(mammaltraits$siteID[i2015])
plots15 <- unique(mammaltraits$plotID[i2015])

# Get local pool for each site

local_pool <- list()

for (site in unique(mam_capture$siteID)) {
  neon_spp <- mam_capture$scientificName[mam_capture$siteID == site]
  regional_spp <- unique(neon_spp)
  regional_spp_id <- unique(mam_capture$taxonID[mam_capture$scientificName %in% regional_spp])
  local_pool <- within(local_pool, assign(site, value=mammaltraits$taxonID %in% regional_spp_id))
}

siteregpoollist_mam15_local <- lapply(local_pool[sites15], function(x) mammaltraits[x, c('weight', 'logweight')])
siteregpoolsp_mam15_local <- lapply(local_pool[sites15], function(x) mammaltraits[x, c('taxonID')])

set.seed(27510)
Ostats_bysite2015 <- Ostats(traits = traits_mam15, plots = factor(mammaltraits$siteID[i2015]), sp = factor(mammaltraits$taxonID[i2015]), reg_pool_traits = siteregpoollist_mam15_local, reg_pool_sp = siteregpoolsp_mam15_local, nperm = N_PERM, stat = STAT, shuffle_weights = FALSE)

o2015 <- ostat2longform(Ostats_bysite2015)

###########################################################################
## Part IV: Regression of overlap statistics on environmental variables. ##
###########################################################################

# Combine overlap statistics and covariates into the same data frame.
o2015 <- left_join(o2015 %>% rename(siteID = site), siteclimatemeans) %>% left_join(richness)

library(betareg)
library(MuMIn)

# Fit full model
reglocal <- betareg(ostat_norm ~ bio1 + chao1 + mpd_z + pc1_productivityheterogeneity + pc2_topographyheterogeneity, link = 'logit', data = o2015 %>% filter(trait == 'logweight'), na.action = 'na.pass')

# Fit all submodels
dredge(reglocal) # Results in best model including only temperature and richness

# Fit best model
reglocal2 <- betareg(ostat_norm ~ bio1 + chao1, link = 'logit', data = o2015 %>% filter(trait == 'logweight'), na.action = 'na.pass')

# Display model fit statistics and parameter estimates.
summary(reglocal2)
confint(reglocal2)

# Check for multicollinearity of predictors.
pred_df <- o2015 %>% filter(trait == 'logweight') %>% select(bio1, chao1, mpd_z, pc1_productivityheterogeneity, pc2_topographyheterogeneity)

library(usdm)
vif(pred_df)
cor(pred_df)

######################################################################################
## Part V: Fit alternative path models relating temperature, richness, and overlap. ##
######################################################################################

library(lavaan)

# logit-transform overlap; z-transform temperature and richness estimator.

semdat <- o2015 %>% 
  filter(trait == 'logweight') %>%
  mutate(logit_ostat = qlogis(ostat_norm),
         bio1_z = (bio1 - mean(bio1, na.rm=T))/sd(bio1, na.rm=T),
         chao1_z = (chao1 - mean(bio1, na.rm=T))/sd(chao1, na.rm=T)) %>%
  dplyr::select(siteID, logit_ostat, bio1_z, chao1_z, mpd_z, pc1_productivityheterogeneity, pc2_topographyheterogeneity)

causal1 <- 
  'logit_ostat ~ bio1_z + chao1_z'

causal2 <-
  'logit_ostat ~ bio1_z + chao1_z
   chao1_z ~ bio1_z'

causal3 <-
  'chao1_z ~ logit_ostat + bio1_z
   logit_ostat ~ bio1_z'

causal4 <-
  'chao1_z ~ logit_ostat
   logit_ostat ~ bio1_z'

causal5 <-
  'logit_ostat ~ chao1_z
   chao1_z ~ bio1_z'

causal6 <-
  'logit_ostat ~~ chao1_z
   chao1_z ~ bio1_z'  

causal7 <-
  'logit_ostat ~~ chao1_z
   logit_ostat ~ bio1_z
   chao1_z ~ bio1_z'  

causal8 <-
  'chao1_z ~ logit_ostat + bio1_z'

causal9 <-
  'logit_ostat ~~ chao1_z
   logit_ostat ~ bio1_z' 

fit1 <- sem(model = causal1, data = semdat)
fit2 <- sem(model = causal2, data = semdat)
fit3 <- sem(model = causal3, data = semdat)
fit4 <- sem(model = causal4, data = semdat)
fit5 <- sem(model = causal5, data = semdat)
fit6 <- sem(model = causal6, data = semdat)
fit7 <- sem(model = causal7, data = semdat)
fit8 <- sem(model = causal8, data = semdat)
fit9 <- sem(model = causal9, data = semdat)

# Extract BIC and R2 from path models.
bics <- rsqs <- rep(NA, 9)

for (i in 1:9) {
  bics[i] <- fitMeasures(get(paste0('fit',i)))['bic']
  rsqs[i] <- max(inspect(get(paste0('fit',i)), 'rsquare'))
}

semfitdat <- data.frame(fit = 1:9, BIC=bics, deltaBIC = bics-min(bics), rsq=rsqs)
