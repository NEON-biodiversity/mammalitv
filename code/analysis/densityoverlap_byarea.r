# Pairwise overlap in absolute area.

source('code/analysis/densityoverlap.r')

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

  

# Run the new pairwise overlap on some of our communities, then compare that to the original.


traits <- traits_mam15
plots <- factor(mam_capture_sitemerge$siteID)
sp <- factor(mam_capture_sitemerge$taxonID)

overlaps_norm_new <- overlaps_norm_old <- matrix(nrow = nlevels(plots), ncol = ncol(traits))


for (s in 1:nlevels(plots)) {
  for (t in 1:ncol(traits)) {
    overlap_norm_st <- try(community_overlap_harmonicwmedian(traits = traits[plots == levels(plots)[s], t], sp = sp[plots == levels(plots)[s]], norm=TRUE), TRUE)
    overlaps_norm_old[s, t] <- if (inherits(overlap_norm_st, 'try-error')) NA else overlap_norm_st
    
    overlap_norm_new_st <- try(community_overlap_absolute(traits = traits[plots == levels(plots)[s], t], sp = sp[plots == levels(plots)[s]], norm=FALSE), TRUE)
    overlaps_norm_new[s, t] <- if (inherits(overlap_norm_new_st, 'try-error')) NA else overlap_norm_new_st
  }
}

comparedat <- data.frame(siteID=unique(plots), ostat_old = overlaps_norm_old[,3], ostat_new = overlaps_norm_new[,3], N = as.numeric(table(plots)))

comparedat <- comparedat %>% left_join(richnessdat) %>% left_join(neonsitedata %>% select(siteID,elevation,bio1,bio4,bio12,bio15,cv_bio1,cv_bio12)) %>% left_join(heterodf) %>% left_join(mammalPDsite)

plot(ostat_new ~ ostat_old, data = comparedat)
plot(ostat_new ~ N, data = comparedat) # This is OK.
plot(log(ostat_new) ~ bio1, data = comparedat)
plot(log(ostat_new) ~ chao1, data = comparedat)
plot(log(ostat_new) ~ mpd_z, data = comparedat)
plot(log(ostat_new) ~ pc1_productivityheterogeneity, data = comparedat)
plot(log(ostat_new) ~ pc2_topographyheterogeneity, data = comparedat)

newlm <- lm(I(log(ostat_new)) ~ bio1 + chao1 + mpd_z + pc1_productivityheterogeneity + pc2_topographyheterogeneity, data = comparedat, na.action = 'na.pass')
summary(newlm)
confint(newlm)

library(MuMIn)
dredge(newlm)
newlm2 <- lm(I(log(ostat_new)) ~ bio1 + chao1 + pc1_productivityheterogeneity + pc2_topographyheterogeneity, data = comparedat, na.action = 'na.pass')
summary(newlm2)
confint(newlm2)
