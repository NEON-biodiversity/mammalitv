source('code/analysis/densityoverlap.r')

x <- c(-12,-10,5,-2,6,10,13)
nindiv <- round(runif(7)*100)
sdnarrow <- c(1.2, 1.3, 0.9, 1.8, 1.3, 1.5, 1)
sdwide <- c(4,3,2.5,6.5,4.5,3.8,3.1)

commnarrow <- data.frame(spid = rep(letters[1:7], times=nindiv),
                         traitval = rnorm(n=sum(nindiv), mean=rep(x, times=nindiv), sd=rep(sdnarrow, times=nindiv)))

commwide <- data.frame(spid = rep(letters[1:7], times=nindiv),
                         traitval = rnorm(n=sum(nindiv), mean=rep(x, times=nindiv), sd=rep(sdwide, times=nindiv)))

pairwise_overlap(a = subset(commnarrow, spid=='a')$traitval, b = subset(commnarrow, spid=='b')$traitval)
pairwise_overlap(a = subset(commwide, spid=='a')$traitval, b = subset(commwide, spid=='b')$traitval)

community_overlap_harmonicwmedian(traits = commnarrow$traitval, sp = commnarrow$spid)
community_overlap_harmonicwmedian(traits = commwide$traitval, sp = commwide$spid)
