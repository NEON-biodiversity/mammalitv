x <- data.frame(effectsizemean=subset(noitv2015goodsites,trait=='logweight')$ostat_norm_localnull_ses,
           effectsizeoverlap=subset(o2015goodsites,trait=='logweight')$ostat_norm_localnull_ses)

plot(x$effectsizemean, -x$effectsizeoverlap)
boxplot(abs(x$effectsizemean)/abs(x$effectsizeoverlap))
abline(0,1)

library(MuMIn)

reglocal <- betareg(ostat_norm ~ bio1 + chao1 + mpd_z + pc1_productivityheterogeneity + pc2_topographyheterogeneity, link = 'logit', data = o2015 %>% filter(trait == 'logweight', !siteID %in% c('DSNY','DELA')), na.action = 'na.pass')

reglocalnoitv <- lm(ostat_norm ~ bio1 + chao1 + mpd_z + ruggedness, data = noitv2015 %>% filter(trait == 'logweight', !siteID %in% c('DSNY','DELA')), na.action = 'na.pass') 

d <- dredge(reglocal)
dnoitv <- dredge(reglocalnoitv)

davg <- model.avg(d)
model.sel(d)


impd <- importance(subset(d, delta < 5))
impdnoitv <- importance(subset(dnoitv, delta < 5))

df_imp <- data.frame(importance = c(impd[1:5], impdnoitv[1:5]), stat = rep(c('Overlap (with ITV)','Mean difference (no ITV)'), each=5),
                     predictor = c(names(impd), names(impdnoitv)))

ggplot(df_imp, aes(x=predictor,y=importance,group=stat, fill=stat)) + geom_bar(stat='identity',position='dodge') +
  theme_john + scale_y_continuous(expand = c(0,0), limits=c(0,1)) +
  scale_x_discrete(labels = c('Temp','Richness','Phy. Dist.','Env. Heterog.','Topo. Heterog.')) +
  theme(legend.position = 'bottom')

ggsave('figs/msfigs/importancebarplot.png', height=6, width=6, dpi=400)
