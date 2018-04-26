# New figures for revised manuscript.
rm(list=ls())
data_path <- '/mnt/research/neon'
source('code/vis/newloadplotdat.r')


# richness ~ temperature

prichtemp <- ggplot(o2015 %>% filter(trait=='logweight'), aes(x=bio1)) + 
  #stat_function(geom='line', fun = fx, args=list(b0 = tempco[1], b1 = tempco[2]), color = 'black', size = 0.8, n=9999) +
  geom_point(aes(y = chao1, color=local_significant), size = 3) +
  geom_text(data=data.frame(bio1=Inf, chao1=Inf, letter='a'), aes(y=chao1,label=letter), size=10, hjust=1, vjust=1) +
  labs(y = 'Species richness (Chao1)', x = parse(text = bioclimnames[1])) +
  theme_john + theme(legend.position = 'none') + csc +
  scale_x_continuous(expand = c(0,0), breaks = c(0,10,20), labels=c(0,10,20), limits=c(-0.5,21.5))

# overlap ~ temperature

# richness ~ overlap

##########
# old

fx <- function(x, b0, b1) 1/(1 + exp(-(b0 + b1 * x)))
fx_loglog <- function(x, b0, b1) exp(-exp(-(b0 + b1 * x)))
tempco <- summary(reglocalbio)$coeff$mean[,1]
chaoco <- summary(reglocalchao)$coeff$mean[,1]
csc <- scale_color_manual(values = c('gray75','black'))

porawtemp <- ggplot(o2015 %>% filter(trait=='logweight'), aes(x=bio1)) + 
  stat_function(geom='line', fun = fx, args=list(b0 = tempco[1], b1 = tempco[2]), color = 'black', size = 0.8, n=9999) +
  geom_point(aes(y = ostat_norm, color=local_significant), size = 3) +
  geom_text(data=data.frame(bio1=Inf, ostat_norm=Inf, letter='b'), aes(y=ostat_norm,label=letter), size=10, hjust=1, vjust=1) +
  labs(y = 'Overlap', x = parse(text = bioclimnames[1])) +
  theme_john + theme(legend.position = 'none') + csc +
  scale_x_continuous(expand = c(0,0), breaks = c(0,10,20), labels=c(0,10,20), limits=c(-0.5,21.5))

porawchao <- ggplot(o2015 %>% filter(trait=='logweight'), aes(x=chao1)) + 
  stat_function(geom='line', fun = fx, args=list(b0 = chaoco[1], b1 = chaoco[2]), color = 'black', size = 0.8, n=9999) +
  geom_point(aes(y = ostat_norm, color=local_significant), size = 3) +
  geom_text(data=data.frame(chao1=Inf, ostat_norm=Inf, letter='c'), aes(y=ostat_norm,label=letter), size=10, hjust=1, vjust=1) +
  labs(y = 'Overlap', x = 'Species richness (Chao1)') +
  theme_john + theme(legend.position = 'none') + csc +
  scale_x_continuous(expand = c(0,0), breaks = c(5,10,15), labels=c(5,10,15), limits=c(2,16))

# regression of chao1~overlap
chaobyoverlm <- lm(chao1 ~ I(qlogis(ostat_norm)), data = o2015 %>% filter(trait=='logweight'))
coef2 <- chaobyoverlm$coefficients
fx2 <- function(x,b0,b1) b0 + b1 * qlogis(x)

pchaobyover <- ggplot(o2015 %>% filter(trait=='logweight'), aes(x=ostat_norm)) + 
  stat_function(geom='line', fun = fx2, args=list(b0 = coef2[1], b1 = coef2[2]), color = 'black', size = 0.8, n=9999) +
  geom_point(aes(y = chao1, color=local_significant), size = 3) +
  geom_text(data=data.frame(chao1=Inf, ostat_norm=Inf, letter='c'), aes(y=chao1,label=letter), size=10, hjust=1, vjust=1) +
  labs(x = 'Overlap', y = 'Species richness (Chao1)') +
  theme_john + theme(legend.position = 'none') + csc +
  scale_x_continuous(expand = c(0,0), breaks = c(0, 0.25, 0.5, 0.75), labels=c(0, 0.25, 0.5, 0.75), limits=c(0,0.91))

library(gridExtra)

png('C:/Users/Q/google_drive/NEON_EAGER/Figures/revisionfigs/threepanelfig_version1.png', height=12, width=5, res=400, units='in')
grid.arrange(prichtemp, porawtemp, porawchao, nrow=3)
dev.off()

png('C:/Users/Q/google_drive/NEON_EAGER/Figures/revisionfigs/threepanelfig_version2.png', height=12, width=5, res=400, units='in')
grid.arrange(prichtemp, porawtemp, pchaobyover, nrow=3)
dev.off()

###########################
# Added 05 April. Scatterplots by productivity (productivity --> overlap as well)

modisdat <- read.csv(file.path(data_path, 'external_data/final_external_data/NEON_modisyearly_withcv.csv'), stringsAsFactors = FALSE)
modis_mean <- modisdat %>% group_by(siteID) %>% summarize_at(-(1:3), mean, na.rm=TRUE)

o2015 <- o2015 %>%
  left_join(modis_mean)

nppexpr <- "paste(\"Net primary productivity (kg m\"^-2,\" y\"^-1,\")\")"
nppreg <- betareg(ostat_norm ~ I(log10(NPP)), link = 'logit', data = o2015 %>% filter(trait=='logweight'))
latco <- latreg$coefficients$mean

porawnpp <- ggplot(o2015 %>% filter(trait=='logweight'), aes(x=NPP)) + 
  #stat_function(geom='line', fun = fx, args=list(b0 = tempco[1], b1 = tempco[2]), color = 'black', size = 0.8, n=9999) +
  geom_point(aes(y = ostat_norm, color=local_significant), size = 3) +
  geom_text(data=data.frame(NPP=Inf, ostat_norm=Inf, letter='d'), aes(y=ostat_norm,label=letter), size=10, hjust=1, vjust=1) +
  labs(y = 'Overlap', x = parse(text=nppexpr)) +
  theme_john + theme(legend.position = 'none') + csc + scale_x_log10(breaks = c(0.2, 0.5, 1, 2, 3))
  #scale_x_continuous(expand = c(0,0), breaks = c(0,10,20), labels=c(0,10,20), limits=c(-0.5,21.5))

nppregchao <- lm(chao1 ~ I(log10(NPP)), data = o2015 %>% filter(trait=='logweight'))
nppco <- nppregchao$coefficients
logline <- function(x, b0, b1) b0 + b1 * log10(x)

prichnpp <- ggplot(o2015 %>% filter(trait=='logweight'), aes(x=NPP)) + 
  #geom_abline(aes(y=chao1), slope = nppco[1], intercept = nppco[2], color = 'black', size = 0.8) +
  stat_function(geom='line', fun = logline, args=list(b0 = nppco[1], b1 = nppco[2]), color = 'black', size = 0.8, n=9999) +
  geom_point(aes(y = chao1, color=local_significant), size = 3) +
  geom_text(data=data.frame(NPP=Inf, chao1=Inf, letter='e'), aes(y=chao1,label=letter), size=10, hjust=1, vjust=1) +
  labs(y = 'Species richness (Chao1)', x = parse(text=nppexpr)) +
  theme_john + theme(legend.position = 'none') + csc + scale_x_log10(breaks = c(0.2, 0.5, 1, 2, 3))
#scale_x_continuous(expand = c(0,0), breaks = c(0,10,20), labels=c(0,10,20), limits=c(-0.5,21.5))

ggsave('C:/Users/Q/google_drive/NEON_EAGER/Figures/revisionfigs/npp_overlap_regression.png', porawnpp, height=6, width=6, dpi=400)
ggsave('C:/Users/Q/google_drive/NEON_EAGER/Figures/revisionfigs/npp_richness_regression.png', prichnpp, height=6, width=6, dpi=400)


# Edit 06 April: create 5 panels separately; SEM will be panel a.
panelb <- ggplot(o2015 %>% filter(trait=='logweight'), aes(x=bio1)) + 
  #stat_function(geom='line', fun = fx, args=list(b0 = tempco[1], b1 = tempco[2]), color = 'black', size = 0.8, n=9999) +
  geom_point(aes(y = chao1, color=local_significant), size = 3) +
  geom_text(data=data.frame(bio1=Inf, chao1=Inf, letter='b'), aes(y=chao1,label=letter), size=10, hjust=1.1, vjust=1.1) +
  labs(y = 'Species richness (Chao1)', x = parse(text = bioclimnames[1])) +
  theme_john + theme(legend.position = 'none') + csc +
  scale_x_continuous(expand = c(0,0), breaks = c(0,10,20), labels=c(0,10,20), limits=c(-0.5,21.5))

# Fix bivariate regression so it's no longer the beta regression. Doesn't matter since it is just for illustrating anyway.
tempreg <- lm(I(qlogis(ostat_norm)) ~ bio1, data = o2015 %>% filter(trait=='logweight'))
tempco2 <- tempreg$coefficients
logitline <- function(b0, b1, x) plogis(b0 + b1 * x)

paneld <- ggplot(o2015 %>% filter(trait=='logweight'), aes(x=bio1)) + 
  #stat_function(geom='line', fun = fx, args=list(b0 = tempco[1], b1 = tempco[2]), color = 'black', size = 0.8, n=9999) +
  stat_function(geom='line', fun = logitline, args=list(b0 = tempco2[1], b1 = tempco2[2]), color = 'black', size = 0.8, n=9999) +
  geom_point(aes(y = ostat_norm, color=local_significant), size = 3) +
  geom_text(data=data.frame(bio1=Inf, ostat_norm=Inf, letter='d'), aes(y=ostat_norm,label=letter), size=10, hjust=1.1, vjust=1.1) +
  labs(y = 'Overlap', x = parse(text = bioclimnames[1])) +
  theme_john + theme(legend.position = 'none') + csc +
  scale_x_continuous(expand = c(0,0), breaks = c(0,10,20), labels=c(0,10,20), limits=c(-0.5,21.5))

panelf <- ggplot(o2015 %>% filter(trait=='logweight'), aes(x=ostat_norm)) + 
  stat_function(geom='line', fun = fx2, args=list(b0 = coef2[1], b1 = coef2[2]), color = 'black', size = 0.8, n=9999) +
  geom_point(aes(y = chao1, color=local_significant), size = 3) +
  geom_text(data=data.frame(chao1=Inf, ostat_norm=Inf, letter='f'), aes(y=chao1,label=letter), size=10, hjust=1.1, vjust=1.1) +
  labs(x = 'Overlap', y = 'Species richness (Chao1)') +
  theme_john + theme(legend.position = 'none') + csc +
  scale_x_continuous(expand = c(0,0), breaks = c(0, 0.25, 0.5, 0.75), labels=c(0, 0.25, 0.5, 0.75), limits=c(0,0.91))

panele <- ggplot(o2015 %>% filter(trait=='logweight'), aes(x=NPP)) + 
  #stat_function(geom='line', fun = fx, args=list(b0 = tempco[1], b1 = tempco[2]), color = 'black', size = 0.8, n=9999) +
  geom_point(aes(y = ostat_norm, color=local_significant), size = 3) +
  geom_text(data=data.frame(NPP=Inf, ostat_norm=Inf, letter='e'), aes(y=ostat_norm,label=letter), size=10, hjust=1.1, vjust=1.1) +
  labs(y = 'Overlap', x = parse(text=nppexpr)) +
  theme_john + theme(legend.position = 'none') + csc + scale_x_log10(breaks = c(0.2, 0.5, 1, 2, 3), limits=c(0.1,3.5), expand=c(0,0))
#scale_x_continuous(expand = c(0,0), breaks = c(0,10,20), labels=c(0,10,20), limits=c(-0.5,21.5))

panelc <- ggplot(o2015 %>% filter(trait=='logweight'), aes(x=NPP)) + 
  #geom_abline(aes(y=chao1), slope = nppco[1], intercept = nppco[2], color = 'black', size = 0.8) +
  stat_function(geom='line', fun = logline, args=list(b0 = nppco[1], b1 = nppco[2]), color = 'black', size = 0.8, n=9999, xlim=c(log10(0.1), log10(4))) +
  geom_point(aes(y = chao1, color=local_significant), size = 3) +
  geom_text(data=data.frame(NPP=Inf, chao1=Inf, letter='c'), aes(y=chao1,label=letter), size=10, hjust=1.1, vjust=1.1) +
  labs(y = 'Species richness (Chao1)', x = parse(text=nppexpr)) +
  theme_john + theme(legend.position = 'none') + csc + scale_x_log10(breaks = c(0.2, 0.5, 1, 2, 3), limits=c(0.1,3.5), expand=c(0,0))

fp <- 'C:/Users/Q/google_drive/NEON_EAGER/Figures/revisionfigs'
ggsave(file.path(fp, 'panelb.png'), panelb, height=5, width=5, dpi=400)
ggsave(file.path(fp, 'panelc.png'), panelc, height=5, width=5, dpi=400)
ggsave(file.path(fp, 'paneld.png'), paneld, height=5, width=5, dpi=400)
ggsave(file.path(fp, 'panele.png'), panele, height=5, width=5, dpi=400)
ggsave(file.path(fp, 'panelf.png'), panelf, height=5, width=5, dpi=400)

# Change letters: 07 April

chaobyoverlm <- lm(chao1 ~ I(qlogis(ostat_norm)), data = o2015 %>% filter(trait=='logweight'))
coef2 <- chaobyoverlm$coefficients
fx2 <- function(x,b0,b1) b0 + b1 * qlogis(x)
tempreg <- lm(I(qlogis(ostat_norm)) ~ bio6, data = o2015 %>% filter(trait=='logweight'))
tempco2 <- tempreg$coefficients
logitline <- function(b0, b1, x) plogis(b0 + b1 * x)

# Also correct for bio6

tl <- "paste(\"Minimum temperature of coldest month (\",degree,\"C)\")" 
tj2 <- theme_john + theme(axis.title=element_text(size=16))

yscrich <- scale_y_continuous(breaks=c(5,10,15), labels=c(5,10,15), limits=c(3,16))

panelb <- ggplot(o2015 %>% filter(trait=='logweight'), aes(x=bio6)) + 
  #stat_function(geom='line', fun = fx, args=list(b0 = tempco[1], b1 = tempco[2]), color = 'black', size = 0.8, n=9999) +
  geom_point(aes(y = chao1), size = 3) +
  geom_text(data=data.frame(bio6=Inf, chao1=Inf, letter='b'), aes(y=chao1,label=letter), size=10, hjust=1.1, vjust=1.1) +
  labs(y = 'Species richness', x = parse(text = tl)) +
  tj2 + theme(legend.position = 'none') + csc +
  yscrich +
  scale_x_continuous(expand = c(0,0), breaks = c(-15,-10,-5,0,5), labels=c(-15,-10,-5,0,5), limits=c(-18.5,8))

panelc <- ggplot(o2015 %>% filter(trait=='logweight'), aes(x=bio6)) + 
  #stat_function(geom='line', fun = fx, args=list(b0 = tempco[1], b1 = tempco[2]), color = 'black', size = 0.8, n=9999) +
  stat_function(geom='line', fun = logitline, args=list(b0 = tempco2[1], b1 = tempco2[2]), color = 'black', size = 0.8, n=9999) +
  geom_point(aes(y = ostat_norm), size = 3) +
  geom_text(data=data.frame(bio6=Inf, ostat_norm=Inf, letter='c'), aes(y=ostat_norm,label=letter), size=10, hjust=1.1, vjust=1.1) +
  labs(y = 'Overlap', x = parse(text = tl)) +
  tj2 + theme(legend.position = 'none') + csc +
  scale_x_continuous(expand = c(0,0), breaks = c(-15,-10,-5,0,5), labels=c(-15,-10,-5,0,5), limits=c(-18.5,8))

paneld <- ggplot(o2015 %>% filter(trait=='logweight'), aes(x=ostat_norm)) + 
  stat_function(geom='line', fun = fx2, args=list(b0 = coef2[1], b1 = coef2[2]), color = 'black', size = 0.8, n=9999) +
  geom_point(aes(y = chao1), size = 3) +
  geom_text(data=data.frame(chao1=Inf, ostat_norm=Inf, letter='d'), aes(y=chao1,label=letter), size=10, hjust=1.1, vjust=1.1) +
  labs(x = 'Overlap', y = 'Species richness') +
  tj2 + theme(legend.position = 'none') + csc +
  yscrich +
  scale_x_continuous(expand = c(0,0), breaks = c(0, 0.25, 0.5, 0.75), labels=c(0, 0.25, 0.5, 0.75), limits=c(0,0.91))

fp <- 'C:/Users/Q/google_drive/NEON_EAGER/Figures/revisionfigs'
ggsave(file.path(fp, '3panelb.png'), panelb, height=5, width=5, dpi=400)
ggsave(file.path(fp, '3panelc.png'), panelc, height=5, width=5, dpi=400)
ggsave(file.path(fp, '3paneld.png'), paneld, height=5, width=5, dpi=400)
