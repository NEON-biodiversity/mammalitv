# Plot and map of overlap by latitude to see if the pattern is striking

rm(list=ls())
source('code/vis/newloadplotdat.r')

# Regression plot

latreg <- betareg(ostat_norm ~ decimalLatitude, link = 'logit', data = o2015 %>% filter(trait=='logweight'))
latco <- latreg$coefficients$mean

porawlat <- ggplot(o2015 %>% filter(trait=='logweight'), aes(x=decimalLatitude)) + 
  stat_function(geom='line', fun = fx, args=list(b0 = latco[1], b1 = latco[2]), color = 'black', size = 0.8, n=9999) +
  geom_point(aes(y = ostat_norm, color=local_significant), size = 3) +
  #geom_text(data=data.frame(bio1=Inf, ostat_norm=Inf, letter='b'), aes(y=ostat_norm,label=letter), size=10, hjust=1, vjust=1) +
  labs(y = 'Overlap', x = 'Latitude') +
  theme_john + theme(legend.position = 'none') + csc +
  scale_x_continuous(expand = c(0,0), breaks = c(30,35,40,45), labels=c(30,35,40,45), limits=c(28,48))

ggsave('C:/Users/Q/Google Drive/NEON_EAGER/Figures/revisionfigs/latitude_regression.png', porawlat, height=6, width=6, dpi=400)

latregchao <- lm(chao1 ~ decimalLatitude, data = o2015 %>% filter(trait=='logweight'))

pchaolat <- ggplot(o2015 %>% filter(trait=='logweight'), aes(x=decimalLatitude)) + 
  geom_point(aes(y = chao1, color=local_significant), size = 3) +
  #geom_text(data=data.frame(bio1=Inf, ostat_norm=Inf, letter='b'), aes(y=ostat_norm,label=letter), size=10, hjust=1, vjust=1) +
  labs(y = 'Species Richness (Chao1)', x = 'Latitude') +
  theme_john + theme(legend.position = 'none') + csc +
  scale_x_continuous(expand = c(0,0), breaks = c(30,35,40,45), labels=c(30,35,40,45), limits=c(28,48))

ggsave('C:/Users/Q/Google Drive/NEON_EAGER/Figures/revisionfigs/latitude_richness_regression.png', pchaolat, height=6, width=6, dpi=400)


# Map

library(ggplot2)

datamap <- ggplot(subset(o2015, decimalLatitude>20 & decimalLatitude<50 & trait == 'logweight'), aes(x=decimalLongitude, y=decimalLatitude)) +
  borders('state', fill='beige') +
  coord_map() + 
  scale_x_continuous(limits = c(-125,-66.2), expand=c(0,0), breaks = c(-120, -105, -90, -75), labels = c('120° W', '105° W', '90° W', '75° W')) +
  scale_y_continuous(limits=c(24.9,50), expand=c(0,0), breaks = c(25, 35, 45), labels = c('25° N', '35° N', '45° N')) +
  theme(panel.border = element_rect(color='black', fill=NA), panel.background=element_blank(), panel.grid=element_blank(), axis.title = element_blank())



mamcs <- scale_fill_gradientn(name = 'Overlap', colours = rev(RColorBrewer::brewer.pal(9, 'RdYlBu')), breaks=c(0,.25,.5,.75,1))

us_map_ovl <- datamap + geom_point(size=5, shape=21, aes(fill = ostat_norm), data=subset(o2015, decimalLatitude>20 & decimalLatitude<50 & trait == 'logweight')) +
  mamcs +
  theme(legend.position=c(0.01,0.01), legend.justification=c(0,0), legend.direction='horizontal')

ggsave('C:/Users/Q/Google Drive/NEON_EAGER/Figures/revisionfigs/mapbysite_ostat.png', us_map_ovl, height=6, width=9, dpi=400)

#usmap(filepath = 'C:/Users/Q/Google Drive/NEON_EAGER/Figures/revisionfigs/mapbysite_ostat.png', fill_var = 'ostat_norm', fill_name = 'Overlap', color_all = rev(RColorBrewer::brewer.pal(9, 'RdYlBu')), puertorico = FALSE, alaska = FALSE, inset_titles = FALSE, parse_name = FALSE, legend_bar_extend = 0.5, font = 'Helvetica', mapdat = subset(o2015, decimalLatitude>20 & decimalLatitude<50 & trait == 'logweight'))
