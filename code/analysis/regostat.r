# Regional O-stats test
# 2015 mammals

# Last modified: 05 Oct 2016

# Still going back and forth between the normalized and non-normalized indices.

traits_bysite <- split(traits_mam15, factor(mam_capture_sitemerge$siteID[i2015]))

plots <- factor(mam_capture_sitemerge$siteID[i2015])

source('code/analysis/densityoverlap.r')

regoverlaps <- nindivlocal <- numeric(length(traits_bysite))


for (i in 1:length(traits_bysite)) {
  tr_i <- traits_mam15[plots == names(siteregpoollist_mam15_iucn)[i], ]
  nindivlocal[i] <- nrow(tr_i)
  regoverlaps[i] <- pairwise_overlap(a = tr_i$logweight, b = siteregpoollist_mam15_iucn[[i]]$logweight, norm=TRUE)[1]
}

# How many indivs in the local site relative to the number of indivs in the local pool
propofregpool <- nindivlocal/sapply(siteregpoollist_mam15_iucn,nrow)

regoverlapdf <- data.frame(siteID=names(siteregpoollist_mam15_iucn), regoverlap=regoverlaps, propofregpool = propofregpool)
regoverlapdf <- transform(regoverlapdf, regoverlap_corrected = regoverlap/propofregpool)

# Test with null model

regoverlap_nulls <- matrix(nrow = 999, ncol=length(traits_bysite))
for (rep in 1:999) {
  for (i in 1:length(traits_bysite)) {
    #tr_i <- traits_mam15[plots == names(siteregpoollist_mam15_iucn)[i], ]
    #nindivlocal[i] <- nrow(tr_i)
    trnull_irep <- sample(x=siteregpoollist_mam15_iucn[[i]]$logweight, size=nindivlocal[i], replace=FALSE)
    regoverlap_nulls[rep,i] <- pairwise_overlap(a = trnull_irep, b = siteregpoollist_mam15_iucn[[i]]$logweight, norm=TRUE)[1]
    
  }
  print(rep)
}

# Null model 2: within-species null model. Only pull the same species from the regional pool that are there.

regoverlap_nulls2 <- matrix(nrow = 999, ncol=length(traits_bysite))
for (rep in 1:999) {
  for (i in 1:length(traits_bysite)) {
    #tr_i <- traits_mam15[plots == names(siteregpoollist_mam15_iucn)[i], ]
    #nindivlocal[i] <- nrow(tr_i)
    sp_i <- factor(mam_capture_sitemerge$taxonID[i2015])[plots == names(siteregpoollist_mam15_iucn)[i]]
    trnull_irep <- c()
    spnames_i <- unique(sp_i)
    sptable_i <- table(sp_i)
    
    for (s in 1:length(spnames_i)) {
      z <- siteregpoollist_mam15_iucn[[i]]$logweight[siteregpoolsp_mam15_iucn[[i]] == as.character(spnames_i[s])]
      if (any(!is.na(z))) trnull_irep <- c(trnull_irep, sample(x=z, size=sum(sp_i == spnames_i[s]), replace=FALSE))
    }  
      #sample(x=siteregpoollist_mam15_iucn[[i]]$logweight, size=nindivlocal[i], replace=FALSE)
    regoverlap_nulls2[rep,i] <- pairwise_overlap(a = trnull_irep, b = siteregpoollist_mam15_iucn[[i]]$logweight, norm=TRUE)[1]

  }
  print(rep)
}




regoverlap_nullmeans <- apply(regoverlap_nulls,2,mean)
regoverlap_nullsds <- apply(regoverlap_nulls,2,sd)
regoverlap_nulllow <- apply(regoverlap_nulls,2,quantile,prob=.025)
regoverlap_nullhigh <- apply(regoverlap_nulls,2,quantile,prob=.975)

regoverlap_ses <- (regoverlaps - regoverlap_nullmeans)/regoverlap_nullsds

regoverlap_nullmeans_nm2 <- apply(regoverlap_nulls2,2,mean)
regoverlap_nullsds_nm2 <- apply(regoverlap_nulls2,2,sd)
regoverlap_nulllow_nm2 <- apply(regoverlap_nulls2,2,quantile,prob=.025)
regoverlap_nullhigh_nm2 <- apply(regoverlap_nulls2,2,quantile,prob=.975)

regoverlap_ses_nm2 <- (regoverlaps - regoverlap_nullmeans_nm2)/regoverlap_nullsds_nm2

regoverlapdf$ses <- regoverlap_ses
regoverlapdf$null_low <- regoverlap_nulllow
regoverlapdf$null_high <- regoverlap_nullhigh

regoverlapdf$ses_nm2 <- regoverlap_ses_nm2
regoverlapdf$null_low_nm2 <- regoverlap_nulllow_nm2
regoverlapdf$null_high_nm2 <- regoverlap_nullhigh_nm2

stattext2 <- rep('neutral',nrow(regoverlapdf))
stattext2[regoverlapdf$regoverlap < regoverlapdf$null_low] <- 'filtered'
stattext2[regoverlapdf$regoverlap > regoverlapdf$null_high] <- 'overdispersed'

stattext2nm2 <- rep('neutral',nrow(regoverlapdf))
stattext2nm2[regoverlapdf$regoverlap < regoverlapdf$null_low_nm2] <- 'filtered'
stattext2nm2[regoverlapdf$regoverlap > regoverlapdf$null_high_nm2] <- 'overdispersed'

regoverlapdf$stattext <- stattext2
regoverlapdf$stattext_nm2 <- stattext2nm2

