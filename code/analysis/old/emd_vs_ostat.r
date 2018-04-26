# Compare pairwise overlap to earth moving distance

x <- mam_capture_sitemerge[i2015, ]
xs <- split(x, x$siteID)

xoverlap <- lapply(xs, function(z) community_overlap_wmedian(traits = z$logweight, sp = z$taxonID))

pairwise_emd <- function(a, b, bw = NULL, n = NULL) {
  
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
  # d <- data.frame(x=da$x, a=da$y, b=db$y)
  
  # If not normalized, multiply each density entry by the length of each vector
  # if (!norm) {
  #   d$a <- d$a * length(a)
  #   d$b <- d$b * length(b)
  # }
  # 
  suppressMessages(require(emdist))
  
  amat <- cbind(da$y, da$x)
  bmat <- cbind(db$y, db$x)
  
  emd(A=amat, B=bmat, dist='euclidean')
  
  # # calculate intersection densities
  # d$w <- pmin(d$a, d$b)
  # 
  # # integrate areas under curves
  # suppressMessages(require(sfsmisc))
  # total <- integrate.xy(d$x, d$a) + integrate.xy(d$x, d$b)
  # intersection <- integrate.xy(d$x, d$w)
  # 
  # # compute overlap coefficient
  # overlap <- 2 * intersection / total
  # overlap_a <- intersection / integrate.xy(d$x, d$a)
  # overlap_b <- intersection / integrate.xy(d$x, d$b)
  # 
  # return(c(overlap = overlap, overlap_a = overlap_a, overlap_b = overlap_b))
  
}

community_emd <- function(traits, sp, norm = TRUE, bw = NULL, n = NULL) {
  sp <- as.character(sp)
  dat <- data.frame(traits=traits, sp=sp, stringsAsFactors = FALSE)
  dat <- dat[complete.cases(dat), ]
  abunds <- table(dat$sp)
  abunds <- abunds[abunds>1]
  dat <- dat[dat$sp %in% names(abunds), ]
  traitlist <- split(dat$traits, dat$sp)
  nspp <- length(traitlist)
  
  if (nspp < 2) return(NA)
  
  eds <- numeric(0)
  abund_pairs <- numeric(0)
  
  for (sp_a in 1:(nspp-1)) {
    for (sp_b in (sp_a+1):nspp) {
      ed <- pairwise_emd(a = traitlist[[sp_a]], b = traitlist[[sp_b]])
      eds <- c(eds, ed)
      abund_pairs <- c(abund_pairs, abunds[sp_a] + abunds[sp_b])
    }
  }
  
  matrixStats::weightedMedian(x = eds, w = abund_pairs)
  
}

xemd <- lapply(xs, function(z) community_emd(traits = z$logweight, sp = z$taxonID))

plot(unlist(xoverlap), unlist(xemd))

# Compare emd to the no-itv case.
xmeandiff <- lapply(xs, function(z) community_overlap_noitv(traits = z$logweight, sp = z$taxonID))

plot(unlist(xmeandiff), unlist(xemd))
abline(0,1,col='red')


plot(unlist(xmeandiff)[-5], 1-unlist(xoverlap)[-5])
logistic <- function(x) 1/(1+exp(-x))
logit <- function(x) log((x/(1-x)))
plot(unlist(xmeandiff)[-5], logit(1 - unlist(xoverlap)[-5]))
abline(0,1,col='red')
