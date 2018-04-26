# function to swap species means around, keeping spreads.

plots <- factor(rep(letters[1:3], each=100))
spmeans1 <- c(1, 1.5, 1.7, 3, 5)
spsds1 <- c(1, 1, 2, 1, 3)
spmeans2 <- c(2, 1, 5, 8, 20)
spsds2 <- c(1.5, 0.5, 2, 1, 3)
sp <- rep(rep(LETTERS[1:5], each=20), 3)
traits1 <- rnorm(300, mean = rep(rep(spmeans1, each=20),3), sd=rep(rep(spsds1, each=20),3))
traits2 <- rnorm(300, mean = rep(rep(spmeans2, each=20),3), sd=rep(rep(spsds2, each=20),3))
traits <- cbind(traits1,traits2)

s <- 1
t <- 1

traits_st <- traits[plots==levels(plots)[s], t]
sp_st <- sp[plots==levels(plots)[s]]

traitmeans <- tapply(traits_st,sp_st,mean)
traitdeviations <- traits_st-traitmeans[sp_st]
#(traitavgdeviation <- tapply(abs(traitdeviations),sp_st,mean))

# Sort the trait means out randomly.
traitmeans_null <- sample(traitmeans)
sp_null <- rep(names(traitmeans_null), table(sp_st))
traits_null <- traitdeviations + traitmeans_null[sp_null]

#tapply(traits_null, names(traits_null), mean)
