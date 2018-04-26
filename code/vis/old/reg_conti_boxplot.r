# Check comparison of local filtering stat between regional and continental filters.
# QDR/NEONITV/2016-11-14/
# Modified 2016-12-02: remove Harvard.
# Modified 2016-12-15: add local

rm(list = ls())
source('code/vis/loadplotdat.r')

# Create jitterplot data
o2015goodsites <- o2015goodsites %>%
  mutate(local_significant = ostat_norm_localnull_ses < ostat_norm_localnull_seslower | ostat_norm_localnull_ses > ostat_norm_localnull_sesupper,
         reg_significant = ostat_norm_regnull_ses < ostat_norm_regnull_seslower | ostat_norm_regnull_ses > ostat_norm_regnull_sesupper)

jitterplotdat <- o2015goodsites %>% 
  filter(trait == 'logweight') %>%
  select(siteID, ostat_norm_localnull_ses, ostat_norm_regnull_ses, local_significant, reg_significant)

jitterplotdat <- with(jitterplotdat, data.frame(siteID=siteID, 
                                                ses = c(ostat_norm_localnull_ses, ostat_norm_regnull_ses),
                                                significant = c(local_significant, reg_significant),
                                                nullmodel = rep(c('Local','Regional'), each=nrow(jitterplotdat))))

jitterplottext <- data.frame(lab = c('More partitioning\nthan expected', 'Neutral', 'More overlap\nthan expected'),
                             x = c(1.5, 1.5, 1.5),
                             y = c(-10, 1, 10))

# Load continental o-stat.
load('C:/Users/Q/Dropbox/neon/code/overlapstat/Ostats_bysite2015_wmedconti.r')

o2015conti <- ostat2longform(Ostats_bysite2015conti)
o2015conti <- o2015conti %>% rename(siteID=site) %>% left_join(neonsitedata) %>% left_join(spstatmeans) %>% left_join(mammalPDsite) %>% left_join(richnessdf) %>% left_join(heterodf) %>% filter(chao1 > 1)
o2015contigoodsites <- filter(o2015conti, !siteID %in% c('DELA','DSNY','HEAL', 'HARV'))

o2015localgoodsites <- filter(orlocal2015, !siteID %in% c('DELA','DSNY','HEAL', 'HARV'))

# compare the jittered plot for continental and regional O-stats.

# Find "significance"

o2015contigoodsites <- o2015contigoodsites %>%
  mutate(local_significant = ostat_norm_localnull_ses < ostat_norm_localnull_seslower | ostat_norm_localnull_ses > ostat_norm_localnull_sesupper,
         reg_significant = ostat_norm_regnull_ses < ostat_norm_regnull_seslower | ostat_norm_regnull_ses > ostat_norm_regnull_sesupper)

jitterplotdatconti <- o2015contigoodsites %>% 
  filter(trait == 'logweight') %>%
  select(siteID, ostat_norm_localnull_ses, ostat_norm_regnull_ses, local_significant, reg_significant)

jitterplotdatconti <- with(jitterplotdatconti, data.frame(siteID=siteID, 
                                                ses = c(ostat_norm_localnull_ses, ostat_norm_regnull_ses),
                                                significant = c(local_significant, reg_significant),
                                                nullmodel = rep(c('Local','Continental'), each=nrow(jitterplotdatconti))))

jitterplotdat <- rbind(jitterplotdat, subset(jitterplotdatconti, nullmodel=='Continental'))

jitterplottext <- data.frame(lab = c('More partitioning\nthan expected', 'Neutral', 'More overlap\nthan expected'),
                             x = c(1.5, 1.5, 1.5),
                             y = c(-10, 1, 10))

pj <- ggplot(jitterplotdat, aes(x=nullmodel,y=ses)) +
  geom_hline(yintercept=0, linetype='dotted', color = 'blue', size=1) +
  geom_jitter(aes(color=significant), height=0, width=0.25) +
  geom_text(aes(x,y,label=lab), data=jitterplottext, family = 'Helvetica') +
  scale_x_discrete(name =  'Null model', labels = c('Local','Regional','Continental')) +
  scale_y_continuous(name = 'Overlap effect size') +
  scale_color_manual(values = c('gray75', 'black')) +
  theme_john + theme(legend.position = c(0.88,0.1))


pjbox <- ggplot(jitterplotdat, aes(x=nullmodel,y=ses)) +
  geom_hline(yintercept=0, linetype='dotted', color = 'blue', size=1) +
  geom_boxplot() +
 # geom_text(aes(x,y,label=lab), data=jitterplottext, family = 'Helvetica') +
  scale_x_discrete(name =  'Null model', labels = c('Local','Regional','Continental')) +
  scale_y_continuous(name = 'Overlap effect size') +
 # scale_color_manual(values = c('gray75', 'black')) +
  theme_john + theme(legend.position = c(0.88,0.1))

ggsave('figs/msfigs/boxplot3nullmodels.png', pjbox, height=5, width=5, dpi=400)


# Medians and quantiles of the effect size distributions
jitterplotdat %>% group_by(nullmodel) %>% summarize(median = median(ses), 
                                                    q025 = quantile(ses, prob=0.025), 
                                                    q975 = quantile(ses, prob=0.975),
                                                    nsigpos = sum(significant & ses>0), 
                                                    nsigneg = sum(significant & ses<0), 
                                                    nnull = sum(!significant))


# Add local

o2015localgoodsites <- o2015localgoodsites %>%
  mutate(local_significant = ostat_norm_localnull_ses < ostat_norm_localnull_seslower | ostat_norm_localnull_ses > ostat_norm_localnull_sesupper,
         reg_significant = ostat_norm_regnull_ses < ostat_norm_regnull_seslower | ostat_norm_regnull_ses > ostat_norm_regnull_sesupper)

jitterplotdatlocal <- o2015localgoodsites %>% 
  filter(trait == 'logweight') %>%
  select(siteID, ostat_norm_localnull_ses, ostat_norm_regnull_ses, local_significant, reg_significant)

jitterplotdatlocal <- with(jitterplotdatlocal, data.frame(siteID=siteID, 
                                                          ses = c(ostat_norm_localnull_ses, ostat_norm_regnull_ses),
                                                          significant = c(local_significant, reg_significant),
                                                          nullmodel = rep(c('Local','Site'), each=nrow(jitterplotdatlocal))))

jitterplotdatall3 <- rbind(subset(jitterplotdatlocal, nullmodel=='Site'), subset(jitterplotdat, nullmodel=='Regional'), subset(jitterplotdatconti, nullmodel=='Continental')) %>% mutate(nullmodel = factor(nullmodel, levels = c('Site','Regional','Continental')))


pjbox <- ggplot(jitterplotdatall3, aes(x=nullmodel,y=ses)) +
  geom_hline(yintercept=0, linetype='dotted', color = 'blue', size=1) +
  geom_boxplot() +
  # geom_text(aes(x,y,label=lab), data=jitterplottext, family = 'Helvetica') +
  scale_x_discrete(name =  'Null model', labels = c('Site','Regional','Continental')) +
  scale_y_continuous(name = 'Overlap effect size') +
  # scale_color_manual(values = c('gray75', 'black')) +
  theme_john + theme(legend.position = c(0.88,0.1))

ggsave('C:/Users/Q/Google Drive/NEON_EAGER/Figures/msfigs/boxplot3nullmodels15Dec.png', pjbox, height=5, width=5, dpi=400)


jitterplotdatall3 %>% group_by(nullmodel) %>% summarize(median = median(ses), 
                                                    q025 = quantile(ses, prob=0.025), 
                                                    q975 = quantile(ses, prob=0.975),
                                                    nsigpos = sum(significant & ses>0), 
                                                    nsigneg = sum(significant & ses<0), 
                                                    nnull = sum(!significant))