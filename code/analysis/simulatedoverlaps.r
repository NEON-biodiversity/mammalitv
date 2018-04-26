# Generate fake data and run overlap stats on them.

sp_means <- c(1, 1, 1)
sp_sds <- c(10, 10, 10)
sp_ns <- c(1000, 1000, 1000)

randomtrait <- rnorm(n = sum(sp_ns), mean = rep(sp_means, sp_ns), sd = rep(sp_sds, sp_ns))
dat <- data.frame(trait = randomtrait, sp = rep(letters[1:length(sp_ns)], sp_ns))

dat <- data.frame(trait=rep(dnorm(seq(-3,3,length.out=100),0,1), times=3), sp=rep(letters[1:3], each=100))

observed_overlap <- community_overlap(traits=dat$trait, sp=dat$sp)

nrep <- 99

null_overlap <- numeric(nrep)

pb <- txtProgressBar(min=0, max=nrep, style=3)

for (rep in 1:nrep) {
  null_overlap[rep] <- community_overlap(traits=dat$trait, sp=sample(dat$sp))
  setTxtProgressBar(pb, rep)
}

close(pb)

quantile(null_overlap, probs= c(0.025, 0.975))
observed_overlap
