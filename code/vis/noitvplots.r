# Plots for mean-only analysis
# 27 Oct 2016

reglocalchaonoharv <- betareg(ostat_norm ~ chao1, link = 'logit', data = o2015 %>% filter(trait == 'logweight', !siteID %in% c('DSNY','DELA','HARV')))

fx <- function(x, b0, b1) 1/(1 + exp(-(b0 + b1 * x)))

porawchao <- ggplot(o2015goodsites %>% filter(trait=='logweight'), aes(x=chao1)) + 
  #geom_segment(aes(y = ostat_norm_localnull_lower, yend = ostat_norm_localnull_upper, xend=chao1), alpha = 0.5, size = 1.5, color = 'skyblue') +
  #geom_segment(aes(y = ostat_norm_regnull_lower, yend = ostat_norm_regnull_upper, xend=chao1), alpha = 0.5, size = 1.5, color = 'goldenrod') +
  geom_point(aes(y = ostat_norm), size = 1.5) +
  stat_function(geom='line', fun = fx, args=list(b0 = 0.556, b1 = -0.1585), color = 'blue', size = 1.5) +
  stat_function(geom='line', fun = fx, args=list(b0 = 1.096, b1 = -0.2478), color = 'red', size = 1.5) +
  labs(y = expression(NO[local]), x = 'Species Richness (Chao1)') +
  theme_john

ggsave('figs/msfigs/regressionharvard.png', height=6, width = 6, dpi=400)


pnoitvchao <- ggplot(noitv2015goodsites %>% filter(trait=='logweight'), aes(x=chao1)) + 
  #geom_segment(aes(y = ostat_norm_localnull_lower, yend = ostat_norm_localnull_upper, xend=chao1), alpha = 0.5, size = 1.5, color = 'skyblue') +
  #geom_segment(aes(y = ostat_norm_regnull_lower, yend = ostat_norm_regnull_upper, xend=chao1), alpha = 0.5, size = 1.5, color = 'goldenrod') +
  geom_point(aes(y = ostat_norm), size = 1.5) +
  #stat_function(geom='line', fun = fx, args=list(b0 = 0.556, b1 = -0.1585), color = 'blue', size = 1.5) +
  #stat_function(geom='line', fun = fx, args=list(b0 = 1.096, b1 = -0.2478), color = 'red', size = 1.5) +
  labs(y = 'Median Pairwise Mean Difference', x = 'Species Richness (Chao1)') +
  theme_john

pnoitvbio <- ggplot(noitv2015goodsites %>% filter(trait=='logweight'), aes(x=bio1)) + 
  #geom_segment(aes(y = ostat_norm_localnull_lower, yend = ostat_norm_localnull_upper, xend=chao1), alpha = 0.5, size = 1.5, color = 'skyblue') +
  #geom_segment(aes(y = ostat_norm_regnull_lower, yend = ostat_norm_regnull_upper, xend=chao1), alpha = 0.5, size = 1.5, color = 'goldenrod') +
  geom_point(aes(y = ostat_norm), size = 1.5) +
  #stat_function(geom='line', fun = fx, args=list(b0 = 0.556, b1 = -0.1585), color = 'blue', size = 1.5) +
  #stat_function(geom='line', fun = fx, args=list(b0 = 1.096, b1 = -0.2478), color = 'red', size = 1.5) +
  labs(y = 'Median Pairwise Mean Difference', x = parse(text=bioclimnames[1])) +
  theme_john


# Jitter plot

library(reshape2)

# Find "significance"

noitv2015goodsites <- noitv2015goodsites %>%
  mutate(local_significant = ostat_norm_localnull_ses < ostat_norm_localnull_seslower | ostat_norm_localnull_ses > ostat_norm_localnull_sesupper,
         reg_significant = ostat_norm_regnull_ses < ostat_norm_regnull_seslower | ostat_norm_regnull_ses > ostat_norm_regnull_sesupper)

jitterplotdat <- noitv2015goodsites %>% 
  filter(trait == 'logweight') %>%
  select(siteID, ostat_norm_localnull_ses, ostat_norm_regnull_ses, local_significant, reg_significant)

jitterplotdat <- with(jitterplotdat, data.frame(siteID=siteID, 
                                                ses = c(ostat_norm_localnull_ses, ostat_norm_regnull_ses),
                                                significant = c(local_significant, reg_significant),
                                                nullmodel = rep(c('Local','Regional'), each=nrow(jitterplotdat))))

jitterplottext <- data.frame(lab = c('More distant\nthan expected', 'Neutral', 'Closer\nthan expected'),
                             x = c(1.5, 1.5, 1.5),
                             y = c(-10, 1, 10))

pj <- ggplot(jitterplotdat, aes(x=nullmodel,y=ses)) +
  geom_hline(yintercept=0, linetype='dotted', color = 'blue', size=1) +
  geom_jitter(aes(color=significant), height=0, width=0.25) +
  geom_text(aes(x,y,label=lab), data=jitterplottext, family = 'Helvetica') +
  scale_x_discrete(name =  'Null model', labels = c('Local','Regional')) +
  scale_y_continuous(name = 'SES of median pairwise distance') +
  scale_color_manual(values = c('gray75', 'black')) +
  theme_john + theme(legend.position = c(0.88,0.1))

ggsave('figs/msfigs/noitv_jitterplot.png', pj, height=5, width=5, dpi=400)
