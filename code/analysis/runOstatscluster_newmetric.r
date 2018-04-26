# Overlap Null Models (without T-stats) run on the mammal data. Run on cluster and split by guild.
# Author: QDR
# Project: NEON ITV
# Created: 11 Oct 2016 *from older script

# 17 Jan: Use Mao-Ning's new metric
# 13 Jan: Update to include new mam_capture data. Not 2016 tho.
# 09 Jan: Remove opportunistic and non-target species. Add option to just do the weight-shuffling null model.
# 20 Dec: Change weighting statistic to harmonic mean of the abundance pairs, rather than geometric mean.
# 14 Dec: Add three levels of region for the null models. Also don't use anything but 2015 and generalists.
# 25 Oct: Added code to get rid of the insectivorous mammals and the individuals not identified to species.
# 13 Oct: Add weighted median.
# 11 Oct: change to run remotely and split by guild. 
# Modified 10 Oct: Add regional O-stats calculation with continent as regional pool.

# Set number of null model iterations here and which stat is used. Also whether the weights should be shuffled.
STAT <- 'harmonic'
N_PERM <- 999
SW <- FALSE 

# Run on cluster: Load data 

# Load neon data
iucn <- read.csv('/mnt/research/neon/external_data/final_external_data/IUCN_mammal_ranges.csv')
load('/mnt/research/neon/final_data/allorganismal2017Jan13.r')
load('/mnt/research/neon/final_data/allsiteplot2016Aug24.r')

mammalTax <- read.csv('/mnt/research/neon/final_data/mammals/NEON_mam_taxonomy.csv')
mammalTraits <- read.csv('/mnt/research/neon/external_data/final_external_data/NEON_miscmammaltraits.csv')

# Source O-stats calc script
source('~/code/ostatcode/Ostats.r')
# Source community density overlap script
source('~/code/ostatcode/densityoverlap.r')

pairwise_overlap <- function(a, b, norm=TRUE, by_area = FALSE, bw = NULL, n = NULL) {
  
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
  
  # Absolute overlap area is just the intersection of the two unnormalized curves.
  if (by_area) return(c(overlap = intersection, overlap_a = intersection, overlap_b = intersection))
  
  # compute overlap coefficient
  overlap <- 2 * intersection / total
  overlap_a <- intersection / integrate.xy(d$x, d$a)
  overlap_b <- intersection / integrate.xy(d$x, d$b)
  
  return(c(overlap = overlap, overlap_a = overlap_a, overlap_b = overlap_b))
  
}

community_overlap_absolute <- function(traits, sp, norm = FALSE, bw = NULL, n = NULL, randomize_weights = FALSE) {
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
      o <- pairwise_overlap(a = traitlist[[sp_a]], b = traitlist[[sp_b]], norm = norm, by_area = TRUE, bw = bw, n = n)
      overlaps <- c(overlaps, o[1])
      harmonic_mean <- 2/(1/abunds[sp_a] + 1/abunds[sp_b])
      abund_pairs <- c(abund_pairs, harmonic_mean)
    }
  }
  
  if (randomize_weights) abund_pairs <- sample(abund_pairs)
  
  matrixStats::weightedMedian(x = overlaps, w = abund_pairs)
  #median(overlaps)
  
}


Ostats <- function(traits, plots, sp, reg_pool_traits, reg_pool_sp, nperm = 99, nullqs = c(0.025, 0.975), shuffle_weights = FALSE) {
  # Required input: a matrix called traits (nrows=n individuals, ncols=n traits), 
  # a vector called plots which is a factor with length equal to nrow(traits),
  # a vector called sp which is a factor with length equal to nrow(traits),
  # a list called reg_pool_traits of which each element is a matrix with ncols=n traits
  # a list called reg_pool_sp of which each element is a factor with length equal to the 
  # number of rows in the corresponding reg_pool_traits matrix.
  
  # Later add some code *here* to check and clean the input, and throw an error
  # if bad input is given. This will not be necessary unless someone besides
  # me ever uses this!
  
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
  #dimnames(regional_overlaps_norm) <- list(levels(plots), dimnames(traits)[[2]])
  #dimnames(regional_overlaps_unnorm) <- list(levels(plots), dimnames(traits)[[2]])
  
  # Calculation of observed O-Stats
  
  print('Calculating observed local O-stats for each community . . .')
  pb <- txtProgressBar(min = 0, max = nlevels(plots), style = 3)
  
  for (s in 1:nlevels(plots)) {
    for (t in 1:ncol(traits)) {
      overlap_norm_st <- try(community_overlap_absolute(traits = traits[plots == levels(plots)[s], t], sp = sp[plots == levels(plots)[s]], norm=FALSE), TRUE)
      overlaps_norm[s, t] <- if (inherits(overlap_norm_st, 'try-error')) NA else overlap_norm_st
      overlap_unnorm_st <- try(community_overlap(traits = traits[plots == levels(plots)[s], t], sp = sp[plots == levels(plots)[s]], norm=FALSE), TRUE)
      overlaps_unnorm[s, t] <- if (inherits(overlap_unnorm_st, 'try-error')) NA else overlap_unnorm_st
    }
    setTxtProgressBar(pb, s)
  }
  
  close(pb)

  # close(pb)
  print('Calculating local and regional null distributions of O-stats . . . ')
  pb <- txtProgressBar(min = 0, max = nperm, style = 3)
  
  # Null model generation and calculation of null O-Stats
  
  # Local null model: generation and calculation done in the same loop
  for (i in 1:nperm) {
    for (s in 1:nlevels(plots)) {
      for (t in 1:ncol(traits)) {
        overlap_norm_sti <- try(community_overlap_absolute(traits = traits[plots == levels(plots)[s], t], sp = sample(sp[plots == levels(plots)[s]]), norm=FALSE), TRUE)
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
        
        regional_overlap_norm_sti <- try(community_overlap_absolute(traits = reg_pool_traits[[s]][null_comm_index, t], sp = reg_pool_sp[[s]][null_comm_index], norm=FALSE), TRUE)
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









# Clean mammal data

library(dplyr)
library(lubridate)
mam_capture <- mutate(mam_capture, 
                      individualandtag = pmin(as.character(individualID), as.character(tagID), na.rm=TRUE),
                      year = year(date)) %>%
  left_join(mammalTax[,c('taxonID','order','taxonProtocolCategory')]) %>% # Added 21 Sep. Gets rid of everything but rodents. 09 Jan: Gets rid of opportunistic and non-target species too.
  filter(order == 'Rodentia', taxonProtocolCategory == 'target', year == 2015)

#mammalGuilds <- mammalTraits %>% 
#  mutate(scientificName = paste(Genus, Species)) %>%
#  dplyr::select(scientificName, Pineda_Main_food)
#
#mammalGuilds <- rbind(mammalGuilds, data.frame(scientificName=c("Dipodomys sp.", "Glaucomys sp.", "Microtus sp.", "Neotoma sp.", 
#                                                                "Perognathus sp.", "Peromyscus sp.", "Reithrodontomys sp."),
#                                               Pineda_Main_food=c('Generalist','Generalist','Herbivore','Herbivore','Generalist','Granivore','Generalist')))


#mam_capture <- left_join(mam_capture, mammalGuilds[-c(89,118,62,63),], by = 'scientificName')

mam_noID <- mam_capture[mam_capture$individualandtag=='', c('year', 'siteID','plotID', 'taxonID','family', 'individualandtag','hindfootLength','earLength','tailLength','totalLength','weight','sex','lifeStage')]

mam_grp <- filter(mam_capture, individualandtag != '') %>%
  group_by(year, siteID, plotID, taxonID, family, individualandtag)

mam_byindiv <- summarize_at(mam_grp, vars(hindfootLength, earLength, tailLength, totalLength, weight), median, na.rm=TRUE)
mam_byindiv_class <- do(mam_grp, data.frame(sex = .$sex[1], lifeStage = names(sort(table(.$lifeStage),decreasing=TRUE))[1]))

mam_byindiv <- cbind(as.data.frame(mam_byindiv), as.data.frame(mam_byindiv_class[,c('sex', 'lifeStage')]))                  
mam_byindiv <-rbind(mam_byindiv, mam_noID)

mam_capture_sitemerge <- merge(mam_byindiv, neonsitedata[,-c(3:5)], sort=FALSE, all.x=TRUE) %>% mutate(logweight = log(weight))

# Here, get rid of the insectivores and the mammals not identified to species.
nottospecies <- c('DPSP','GLSP','MISP','NESP','PESP','PGSP','RESP')
mam_capture_sitemerge <- mam_capture_sitemerge %>% filter(!(taxonID %in% nottospecies) &!(siteID %in% c('DSNY','DELA','HEAL'))) # Gets rid of ~700.

# Split into year+guild combination
#i2015 <- mam_capture_sitemerge$year == 2015

# Convert adult mammal data to trait-only format
traits_mam15 <- mam_capture_sitemerge[, c('hindfootLength', 'weight', 'logweight')]

# Define IUCN Regional Pools

mamnames <- unique(mam_capture$scientificName) 

nomatchspp <- mamnames[!mamnames %in% iucn$binomial]
# Get rid of the individual Tamias alpinus which should not be in regional species pools
# nomatchspp <- nomatchspp[!nomatchspp %in% c('Tamias alpinus')]

sites15 <- unique(mam_capture_sitemerge$siteID)
plots15 <- unique(mam_capture_sitemerge$plotID)

# Get local pool for each site

local_pool <- list()

for (site in unique(mam_capture$siteID)) {
  neon_spp <- mam_capture$scientificName[mam_capture$siteID == site]
  regional_spp <- unique(c(neon_spp, nomatchspp))
  regional_spp_id <- unique(mam_capture$taxonID[mam_capture$scientificName %in% regional_spp])
  local_pool <- within(local_pool, assign(site, value=mam_capture_sitemerge$taxonID %in% regional_spp_id))
}

siteregpoollist_mam15_local <- lapply(local_pool[sites15], function(x) mam_capture_sitemerge[x, c('hindfootLength', 'weight', 'logweight')])
siteregpoolsp_mam15_local <- lapply(local_pool[sites15], function(x) mam_capture_sitemerge[x, c('taxonID')])


# Get IUCN regional pool for each site

reg_pool <- list()

for (site in unique(mam_capture$siteID)) {
  iucn_spp <- as.character(iucn$binomial[iucn[,site]])
  neon_spp <- mam_capture$scientificName[mam_capture$siteID == site]
  regional_spp <- unique(c(iucn_spp, neon_spp, nomatchspp))
  regional_spp_id <- unique(mam_capture$taxonID[mam_capture$scientificName %in% regional_spp])
  reg_pool <- within(reg_pool, assign(site, value=mam_capture_sitemerge$taxonID %in% regional_spp_id))
}

siteregpoollist_mam15_iucn <- lapply(reg_pool[sites15], function(x) mam_capture_sitemerge[x, c('hindfootLength', 'weight', 'logweight')])
siteregpoolsp_mam15_iucn <- lapply(reg_pool[sites15], function(x) mam_capture_sitemerge[x, c('taxonID')])

# Define entire continent as regional pool for each of the local and regional stats.
continentpool15 <- lapply(reg_pool[sites15], function(x) mam_capture_sitemerge[, c('hindfootLength', 'weight', 'logweight')])
continentsp15 <- lapply(reg_pool[sites15], function(x) mam_capture_sitemerge[, c('taxonID')])

# RUN THE OSTATS

task <- as.numeric(Sys.getenv('PBS_ARRAYID'))
print(task)

set.seed(27510)

if (task == 1)
{Ostats_bysite2015local <- Ostats(traits = traits_mam15, plots = factor(mam_capture_sitemerge$siteID), sp = factor(mam_capture_sitemerge$taxonID), reg_pool_traits = siteregpoollist_mam15_local, reg_pool_sp = siteregpoolsp_mam15_local, nperm = N_PERM, stat = STAT, shuffle_weights = SW)
save(Ostats_bysite2015local, file = '~/data/Ostats09Jan/Ostats_bysite2015_harmoniclocal_newmetric.r')}

if (task == 2)
{Ostats_bysite2015iucn <- Ostats(traits = traits_mam15, plots = factor(mam_capture_sitemerge$siteID), sp = factor(mam_capture_sitemerge$taxonID), reg_pool_traits = siteregpoollist_mam15_iucn, reg_pool_sp = siteregpoolsp_mam15_iucn, nperm = N_PERM, stat = STAT, shuffle_weights = SW)
save(Ostats_bysite2015iucn, file = '~/data/Ostats09Jan/Ostats_bysite2015_harmoniciucn_newmetric.r')}

if (task == 3)
{Ostats_bysite2015conti <- Ostats(traits = traits_mam15, plots = factor(mam_capture_sitemerge$siteID), sp = factor(mam_capture_sitemerge$taxonID), reg_pool_traits = continentpool15, reg_pool_sp = continentsp15, nperm = N_PERM, stat = STAT, shuffle_weights = SW)
save(Ostats_bysite2015conti, file = '~/data/Ostats09Jan/Ostats_bysite2015_harmonicconti_newmetric.r')}
