# T-statistics calculation. Tweaked to use weighted means at the community level.
# Author: Adrien Taudiere, edited by QDR
# Project: NEON ITV
# Created: 26 Sep 2016
# Last modified: 28 Sep 2016

# Added 28 Sep 2016: add pairwise overlap stats too.
# Added 29 Sep 2016: add error catching function to return NA if the pairwise overlap can't be calculated for a community.

QTstats <- function (traits, ind.plot, sp, SE = 0, reg.pool = NULL, SE.reg.pool = NULL, 
  nperm = 99, printprogress = TRUE) 
{
  if (sum(is.na(traits)) > 0) {
    warning(paste("This function excludes", sum(is.na(traits)), 
      "Na values", sep = " "))
  }
  if (!is.matrix(traits)) {
    traits <- as.matrix(traits)
  }
  names_sp_ind.plot <- as.factor(paste(sp, ind.plot, sep = "@"))
  Tplosp <- unlist(strsplit(levels(names_sp_ind.plot), split = "@"))[2 * 
    (1:nlevels(names_sp_ind.plot))]
  names(Tplosp) <- levels(names_sp_ind.plot)
  if (!is.null(nperm)) {
    if (nperm == 0) {
      nperm = NULL
    }
  }
  if (length(SE) == 1) {
    SE <- rep(SE, times = ncol(traits))
  }
  mean_IP <- matrix(nrow = nlevels(names_sp_ind.plot), ncol = ncol(traits))
  rownames(mean_IP) = levels(names_sp_ind.plot)
  n_IP <- mean_IP
  mean_PC <- matrix(nrow = nlevels(ind.plot), ncol = ncol(traits))
  var_IP <- matrix(nrow = nlevels(names_sp_ind.plot), ncol = ncol(traits))
  var_PC <- matrix(nrow = nlevels(ind.plot), ncol = ncol(traits))
  var_CR <- vector()
  var_IC <- matrix(nrow = nlevels(ind.plot), ncol = ncol(traits))
  var_PR <- vector()
  var_IR <- vector()
  T_IP.IC <- matrix(nrow = nlevels(ind.plot), ncol = ncol(traits))
  T_IC.IR <- matrix(nrow = nlevels(ind.plot), ncol = ncol(traits))
  T_PC.PR <- matrix(nrow = nlevels(ind.plot), ncol = ncol(traits))
  
  # Added 28 Sep:
  overlaps <- matrix(nrow = nlevels(ind.plot), ncol = ncol(traits))
  
  library(plyr)
  library(Hmisc)
  
  for (t in 1:ncol(traits)) {
    mean_IP[, t] <- tapply(traits[, t], names_sp_ind.plot, 
      mean, na.rm = T) # Mean values for each species plot combination (don't weight)
	  
	n_IP[, t] <- tapply(traits[, t], names_sp_ind.plot, function(x) sum(!is.na(x)))
	
	# The following line was edited by Q, 26 Sep 2016
	dat_PC <- data.frame(mean_IP = mean_IP[, t], n_IP = n_IP[,t], Tplosp = Tplosp)
	
    #mean_PC[, t] <- tapply(mean_IP[, t], Tplosp, mean, na.rm = T) # Mean values for each plot (mean of species means) (WEIGHT)
	
	# The following line was edited by Q, 26 Sep 2016
	mean_PC[, t] <- as.numeric(plyr::ddply(dat_PC, .(Tplosp), function(z) data.frame(wm = wtd.mean(x = z$mean_IP, weights = z$n_IP)))$wm)
	
    var_IP[, t] <- tapply(traits[, t], names_sp_ind.plot, 
      var, na.rm = T) #  Variances for each species plot combination
    
	#var_PC[, t] <- tapply(mean_IP[, t], Tplosp, var, na.rm = T)
	
	# The following line was edited by Q, 26 Sep 2016
	var_PC[, t] <- as.numeric(plyr::ddply(dat_PC, .(Tplosp), function(z) data.frame(wv = wtd.var(x = z$mean_IP, weights = z$n_IP)))$wv)
	
    var_CR[t] <- var(mean_PC[, t], na.rm = T) # Variance of population means (don't weight)
    var_IC[, t] <- tapply(traits[, t], ind.plot, var, na.rm = T) # Variance for each plot, ignoring species (don't weight)
	
    #var_PR[t] <- var(as.vector(mean_IP[, t]), na.rm = T) # Variance of all species means, ignoring abundances (WEIGHT)
	
	# The following line was edited by Q, 26 Sep 2016
	var_PR[t] <- wtd.var(x = as.vector(mean_IP[, t]), weights = as.vector(n_IP[,t]))
	
    var_IR[t] <- var(traits[, t], na.rm = T) # Variance for regional species pool (entire trait data set), ignoring species (don't weight)
    for (s in 1:nlevels(ind.plot)) {
      #T_IP.IC[s, t] <- mean(var_IP[grepl(levels(ind.plot)[s], 
      #  Tplosp), t], na.rm = T)/var_IC[s, t]
	
      # The following line was edited by Q, 26 Sep 2016	
	  T_IP.IC[s, t] <- wtd.mean(x = var_IP[grepl(levels(ind.plot)[s], Tplosp), t], 
								weights = n_IP[grepl(levels(ind.plot)[s], Tplosp), t])/var_IC[s, t]
		
      T_IC.IR[s, t] <- var_IC[s, t]/var_IR[t]
      T_PC.PR[s, t] <- var_PC[s, t]/var_PR[t]
	  
	  # Added 28 Sep:
	  overlap_st <- try(community_overlap(traits = traits[ind.plot == levels(ind.plot)[s], t], sp = sp[ind.plot == levels(ind.plot)[s]]), TRUE)
	  overlaps[s, t] <- if (inherits(overlap_st, 'try-error')) NA else overlap_st
    }
		
  }
  if (is.numeric(nperm)) {
    var_IP_nm1 <- array(dim = c(nperm, ncol(traits), nrow = length(Tplosp)))
    var_PC_nm2sp <- array(dim = c(nperm, ncol(traits), nlevels(ind.plot)))
    var_IC_nm1 <- array(dim = c(nperm, ncol(traits), nlevels(ind.plot)))
    var_IC_nm2 <- array(dim = c(nperm, ncol(traits), nlevels(ind.plot)))
    var_PR_nm2sp <- array(dim = c(nperm, ncol(traits)))
    var_IR_nm2 <- array(dim = c(nperm, ncol(traits)))
    mean_IP_nm2sp <- array(dim = c(nperm, ncol(traits), 
      length(Tplosp)))
	n_IP_nm2sp <- mean_IP_nm2sp  
    mean_PC_nm2sp <- array(dim = c(nperm, ncol(traits), 
      nlevels(ind.plot)))
    traits.nm1 <- list()
    traits.nm2 <- list()
    traits.nm2sp <- list()
    T_IP.IC_nm1 <- array(dim = c(nperm, ncol(traits), nlevels(ind.plot)))
    T_IC.IR_nm2 <- array(dim = c(nperm, ncol(traits), nlevels(ind.plot)))
    T_PC.PR_nm2sp <- array(dim = c(nperm, ncol(traits), 
      nlevels(ind.plot)))
	overlaps_nm <- array(dim = c(nperm, ncol(traits), nlevels(ind.plot)))  
	  
    if (!is.null(reg.pool)) {
      isnullregpool <- FALSE
    }
    else {
      isnullregpool <- TRUE
    }
    if (is.null(reg.pool)) {
      reg.pool <- rep(list(traits), nlevels(ind.plot))
    }
    if (is.data.frame(reg.pool) | is.matrix(reg.pool)) {
      reg.pool <- rep(list(reg.pool), nlevels(ind.plot))
    }
    if (is.null(SE.reg.pool)) {
      SE.reg.pool <- SE
    }
    if (length(SE.reg.pool) == 1) {
      SE.reg.pool <- rep(SE.reg.pool, times = ncol(traits))
    }
    if (is.vector(SE.reg.pool)) {
      SE.reg.pool <- rep(list(SE.reg.pool), nlevels(ind.plot))
    }
    if (length(reg.pool) != nlevels(ind.plot)) {
      stop("reg.pool need to be either a matrix or a list of length equal to the number of communities")
    }
    if (length(SE) != ncol(traits) & length(SE) != 1) {
      stop("The vector SE need to have a length of one or equal to the number of traits")
    }
    if (!isnullregpool) {
      if (length(SE.reg.pool) != length(reg.pool)) {
        stop("The vector SE.reg.pool need to have the same dimension as reg.pool")
      }
    }
    if (printprogress == T) {
      print("creating null models")
    }
    for (t in 1:ncol(traits)) {
      traits.nm1[[t]] <- list()
      for (s in 1:nlevels(ind.plot)) {
        traits.nm1[[t]][[s]] <- list()
        for (i in 1:nperm) {
          trait.intern <- traits[, t]
          if (SE[t] != 0) {
            trait.intern <- rnorm(length(trait.intern), 
              mean = trait.intern, sd = SE[t])
          }
          if (length(traits[ind.plot == levels(ind.plot)[s], 
            t]) != 1) {
            perm_ind.plot1 <- sample(trait.intern[ind.plot == 
              levels(ind.plot)[s]], table(ind.plot)[s])
            traits.nm1[[t]][[s]][[i]] <- perm_ind.plot1
          }
          else {
            traits.nm1[[t]][[s]][[i]] <- "NA"
          }
        }
      }
      if (printprogress == T) {
        print(paste(round(t/ncol(traits)/3 * 100, 2), 
          "%"))
      }
      else {
      }
    }
    for (t in 1:ncol(traits)) {
      traits.nm2[[t]] <- list()
      for (s in 1:nlevels(ind.plot)) {
        traits.nm2[[t]][[s]] <- list()
        for (i in 1:nperm) {
          trait.intern <- reg.pool[[s]][, t]
          SE.reg.pool.intern <- SE.reg.pool[[s]][t]
          if (SE.reg.pool.intern != 0) {
            trait.intern <- rnorm(length(trait.intern), 
              mean = trait.intern, sd = SE.reg.pool.intern)
          }
          perm_ind.plot2 <- sample(trait.intern, table(ind.plot)[s])
          traits.nm2[[t]][[s]][[i]] <- perm_ind.plot2
        }
      }
      if (printprogress == T) {
        print(paste(round(33.3 + t/ncol(traits)/3 * 
          100, 2), "%"))
      }
      else {
      }
    }
    traits_by_sp <- apply(traits, 2, function(x) tapply(x, 
      names_sp_ind.plot, mean))
    traits_by_pop <- traits_by_sp[match(names_sp_ind.plot, 
      rownames(traits_by_sp)), ]
    if (!is.matrix(traits_by_pop)) {
      traits_by_pop <- as.matrix(traits_by_pop)
    }
    for (t in 1:ncol(traits)) {
      traits.nm2sp[[t]] <- list()
      for (s in 1:nlevels(ind.plot)) {
        traits.nm2sp[[t]][[s]] <- list()
        for (i in 1:nperm) {
          perm_ind.plot2sp <- sample(traits_by_pop[, 
            t], table(ind.plot)[s])
          traits.nm2sp[[t]][[s]][[i]] <- perm_ind.plot2sp
        }
      }
      if (printprogress == T) {
        print(paste(round(66.6 + t/ncol(traits)/3 * 
          100, 2), "%"))
      }
      else {
      }
    }
    if (printprogress == T) {
      print("calculation of Tstats using null models")
    }
    yy <- length(names_sp_ind.plot)
    for (t in 1:ncol(traits)) {
      for (i in 1:nperm) {
        mean_IP_nm2sp[i, t, ] <- tapply(unlist(traits.nm2sp[[t]])[(1 + 
          (i - 1) * yy):(i * yy)], names_sp_ind.plot, 
          function(x) mean(x, na.rm = T))
				
		# ADD n_IP line and data frame creation line  
		
		n_IP_nm2sp[i, t, ] <- tapply(unlist(traits.nm2sp[[t]])[(1 + 
          (i - 1) * yy):(i * yy)], names_sp_ind.plot, 
          function(x) sum(!is.na(x)))
		
		dat_PC <- data.frame(mean_IP = mean_IP_nm2sp[i, t, ], n_IP = n_IP_nm2sp[i, t, ], Tplosp = Tplosp)
		
        #mean_PC_nm2sp[i, t, ] <- tapply(mean_IP_nm2sp[i, 
        #  t, ], Tplosp, mean, na.rm = T) # CHANGE THIS
		
		mean_PC_nm2sp[i, t, ] <- as.numeric(plyr::ddply(dat_PC, .(Tplosp), function(z) data.frame(wm = wtd.mean(x = z$mean_IP, weights = z$n_IP)))$wm)
		
      }
      if (printprogress == T) {
        print(paste(round(t/ncol(traits)/3 * 100, 2), 
          "%"))
      }
      else {
      }
    }
    for (t in 1:ncol(traits)) {
      for (i in 1:nperm) {
        var_IP_nm1[i, t, ] <- tapply(unlist(traits.nm1[[t]])[(1 + 
          (i - 1) * yy):(i * yy)], names_sp_ind.plot, 
          function(x) var(x, na.rm = T))
        #var_PC_nm2sp[i, t, ] <- tapply(mean_IP_nm2sp[i, 
        #  t, ], Tplosp, var, na.rm = T) # CHANGE THIS
		
		dat_PC <- data.frame(mean_IP = mean_IP_nm2sp[i, t, ], n_IP = n_IP_nm2sp[i, t, ], Tplosp = Tplosp)
		var_PC_nm2sp[i, t, ] <- as.numeric(plyr::ddply(dat_PC, .(Tplosp), function(z) data.frame(wv = wtd.var(x = z$mean_IP, weights = z$n_IP)))$wv)	
		  
        var_IC_nm1[i, t, ] <- tapply(unlist(traits.nm1[[t]])[(1 + 
          (i - 1) * yy):(i * yy)], ind.plot, function(x) var(x, 
          na.rm = T))
        var_IC_nm2[i, t, ] <- tapply(unlist(traits.nm2[[t]])[(1 + 
          (i - 1) * yy):(i * yy)], ind.plot, function(x) var(x, 
          na.rm = T))
        #var_PR_nm2sp[i, t] <- var(as.vector(mean_IP_nm2sp[i, 
        #  t, ]), na.rm = T) # CHANGE THIS
		 
		var_PR_nm2sp[i, t] <- wtd.var(x = as.vector(mean_IP_nm2sp[i, t, ]), weights = as.vector(n_IP_nm2sp[i, t, ])) 
		  
        var_IR_nm2[i, t] <- var(unlist(traits.nm2[[t]])[(1 + 
          (i - 1) * yy):(i * yy)], na.rm = T)
      }
      if (printprogress == T) {
        print(paste(round(33.3 + t/ncol(traits)/3 * 
          100, 2), "%"))
      }
      else {
      }
    }
    for (t in 1:ncol(traits)) {
      for (i in 1:nperm) {
        for (s in 1:nlevels(ind.plot)) {
          #T_IP.IC_nm1[i, t, s] <- mean(var_IP_nm1[i, 
          #  t, grepl(levels(ind.plot)[s], Tplosp)], 
          #  na.rm = T)/var_IC_nm1[i, t, s] # CHANGE THIS
		  
		  T_IP.IC_nm1[i, t, s] <- wtd.mean(x = var_IP_nm1[i, t, grepl(levels(ind.plot)[s], Tplosp)], 
										   weights = n_IP_nm2sp[i, t, grepl(levels(ind.plot)[s], Tplosp)])/var_IC_nm1[i, t, s]
		  
          T_IC.IR_nm2[i, t, s] <- var_IC_nm2[i, t, s]/var_IR_nm2[i, 
            t]
          T_PC.PR_nm2sp[i, t, s] <- var_PC_nm2sp[i, 
            t, s]/var_PR_nm2sp[i, t]
			
		  overlap_its <- try(community_overlap(traits = unlist(traits.nm1[[t]])[(1 + 
            (i - 1) * yy):(i * yy)][ind.plot == levels(ind.plot)[s]], sp = sp[ind.plot == levels(ind.plot)[s]]) ,TRUE)	
		  overlaps_nm[i, t, s] <- if (inherits(overlap_its, 'try-error')) NA else overlap_its
        }
      }
      if (printprogress == T) {
        print(paste(round(66.6 + t/ncol(traits)/3 * 
          100, 2), "%"))
      }
      else {
      }
    }
  }
  colnames(T_IP.IC) <- colnames(traits)
  colnames(T_IC.IR) <- colnames(traits)
  colnames(T_PC.PR) <- colnames(traits)
  if (is.numeric(nperm)) {
    colnames(T_IP.IC_nm1) <- colnames(traits)
    colnames(T_IC.IR_nm2) <- colnames(traits)
    colnames(T_PC.PR_nm2sp) <- colnames(traits)
  }
  rownames(T_IP.IC) <- levels(as.factor(Tplosp))
  rownames(T_IC.IR) <- levels(as.factor(Tplosp))
  rownames(T_PC.PR) <- levels(as.factor(Tplosp))
  res <- list()
  res$Tstats <- list()
  res$Tstats$T_IP.IC <- T_IP.IC
  res$Tstats$T_IC.IR <- T_IC.IR
  res$Tstats$T_PC.PR <- T_PC.PR
  res$variances <- list()
  res$variances$var_IP <- var_IP
  res$variances$var_PC <- var_PC
  res$variances$var_CR <- var_CR
  res$variances$var_IC <- var_IC
  res$variances$var_PR <- var_PR
  res$variances$var_IR <- var_IR
  if (is.numeric(nperm)) {
    res$variances$var_IP_nm1 <- var_IP_nm1
    res$variances$var_PC_nm2sp <- var_PC_nm2sp
    res$variances$var_IC_nm1 <- var_IC_nm1
    res$variances$var_IC_nm2 <- var_IC_nm2
    res$variances$var_PR_nm2sp <- var_PR_nm2sp
    res$variances$var_IR_nm2 <- var_IR_nm2
    res$Tstats$T_IP.IC_nm <- T_IP.IC_nm1
    res$Tstats$T_IC.IR_nm <- T_IC.IR_nm2
    res$Tstats$T_PC.PR_nm <- T_PC.PR_nm2sp
  }
  else {
  }
  res$traits <- traits
  res$ind.plot <- ind.plot
  res$sp <- sp
  res$sites_richness <- table(ind.plot)
  res$namestraits <- colnames(traits)
  res$call <- match.call()
  res$overlaps <- overlaps
  res$overlaps_nm <- overlaps_nm
  class(res) <- "Tstats"
  invisible(res)
}
