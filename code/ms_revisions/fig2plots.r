# Figure 2 panels a and b
# Revised version created 19 Jun 2017


rm(list=ls())
source('code/vis/newloadplotdat.r')

# Panel a: Map

datamap <- ggplot(subset(neonsitedata, decimalLatitude>20 & decimalLatitude<50), aes(x=decimalLongitude, y=decimalLatitude)) +
  borders('state', fill='beige') +
  coord_map() + 
  scale_x_continuous(limits = c(-125,-66.2), expand=c(0,0), breaks = c(-120, -105, -90, -75), labels = c('120° W', '105° W', '90° W', '75° W')) +
  scale_y_continuous(limits=c(24.9,50), expand=c(0,0), breaks = c(25, 35, 45), labels = c('25° N', '35° N', '45° N')) +
  theme(panel.border = element_rect(color='black', fill=NA), panel.background=element_blank(), panel.grid=element_blank(), axis.title = element_blank())

mamcs <- scale_fill_gradientn(name = 'Overlap', colours = RColorBrewer::brewer.pal(9, 'RdYlBu'), breaks=c(0,.25,.5,.75,1))

us_map_ovl <- datamap + geom_point(size=5, shape=21, aes(fill = ostat_norm), data=subset(o2015, decimalLatitude>20 & decimalLatitude<50 & trait == 'logweight')) +
  mamcs +
  theme(legend.position=c(0.01,0.01), legend.justification=c(0,0), legend.direction='horizontal')

ggsave('C:/Users/Q/google_drive/NEON_EAGER/Figures/revisionfigs/mapbysite_ostat_19Jun.png', us_map_ovl, height=4.15, width=9*5/6, dpi=400)

# Panel b: Regression plot

latreg <- betareg(ostat_norm ~ decimalLatitude, link = 'logit', data = o2015 %>% filter(trait=='logweight'))
latco <- latreg$coefficients$mean

porawlat <- ggplot(o2015 %>% filter(trait=='logweight'), aes(x=decimalLatitude)) + 
  stat_function(geom='line', fun = fx, args=list(b0 = latco[1], b1 = latco[2]), color = 'black', size = 0.8, n=9999) +
  geom_point(aes(y = ostat_norm), color = 'black', size = 3) +
  #geom_text(data=data.frame(bio1=Inf, ostat_norm=Inf, letter='b'), aes(y=ostat_norm,label=letter), size=10, hjust=1, vjust=1) +
  labs(y = 'Overlap', x = 'Latitude') +
  theme_john + theme(legend.position = 'none') + 
  scale_x_continuous(expand = c(0,0), breaks = c(30,35,40,45), labels=paste0(c(30,35,40,45), '° N'), limits=c(28,48))

ggsave('C:/Users/Q/google_drive/NEON_EAGER/Figures/revisionfigs/latitude_regression_19jun.png', porawlat, height=4.15, width=4.15, dpi=400)
