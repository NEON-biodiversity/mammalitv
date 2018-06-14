source('~/GitHub/NEON_repos/mammalitv/code/vis/newloadplotdat.r')

o2015 <- filter(o2015,trait=='logweight')

summary(lm(chao1 ~ bio6, data=o2015))
summary(lm(qlogis(ostat_norm) ~ bio6, data=o2015))
summary(lm(chao1 ~ qlogis(ostat_norm), data=o2015))
