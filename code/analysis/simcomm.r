# Function to simulate data on which to calculate overlap metrics
# So far, can do normal and uniform distributions.

sim_comm <- function(n_spp, abunds = rep(10, n_spp), distributions = rep('norm', n_spp), means = NULL, sds = NULL, a = NULL, b = NULL) {
  trait_values <- list()
  for (i in 1:n_spp) {
    fx <- get(paste0('r',distributions[i]))
    if (distributions[i] == 'norm')
      trait_values[[i]] <- fx(n = abunds[i], mean = means[i], sd = sds[i])
    if (distributions[i] == 'unif')
      trait_values[[i]] <- fx(n = abunds[i], min = a[i], max = b[i])
  }
  data.frame(species = rep(LETTERS[1:n_spp], times = abunds),
             trait = unlist(trait_values))
}

# Calculation of old and new metrics on different simulated datasets.
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


comm1 <- sim_comm(n_spp = 3, abunds = c(10,10,100), means = c(1,2,1), sds = c(1,1,1))
community_overlap(comm1$trait, comm1$species, norm = TRUE)
community_overlap_absolute(comm1$trait, comm1$species, norm = FALSE)

# Test what happens when the abundance of species 3 becomes greater and greater.

ns <- 3:1000

oldmetrics <- newmetrics <- numeric(length(ns))

for (n in 1:length(ns)) {
  comm_n <- sim_comm(n_spp = 3, abunds = c(10,10,ns[n]), means = c(1,2,1), sds = c(1,1,1))
  oldmetrics[n] <- community_overlap(comm_n$trait, comm_n$species, norm = TRUE)
  newmetrics[n] <- community_overlap_absolute(comm_n$trait, comm_n$species, norm = FALSE)
}

plot(ns, oldmetrics)
plot(ns, newmetrics)

# Result: As one species takes over the community in abundance, the overlap begins to max out at the total number of other individuals.