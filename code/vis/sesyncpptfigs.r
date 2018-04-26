# Scatterplots with white on black theme for PPTs.

bktheme <- theme(text = element_text(color = 'white'),
                 panel.grid = element_blank(),
                 panel.border = element_rect(color = 'gray50', fill = 'transparent'),
                 panel.background = element_rect(fill = 'black'),
                 plot.background = element_rect(fill = 'black', color = 'black'),
                 axis.text.x = element_text(size = 12, color = 'white'), 
                 axis.title = element_text(size = 16, color = 'white'), 
                 axis.text.y = element_blank(),
                 axis.ticks.y = element_blank(),
                 legend.position = 'none',
                 strip.background = element_blank(),
                 strip.text = element_text(family = 'Helvetica', color = 'white'))


panelb <- ggplot(o2015 %>% filter(trait=='logweight'), aes(x=bio6)) + 
  #stat_function(geom='line', fun = fx, args=list(b0 = tempco[1], b1 = tempco[2]), color = 'black', size = 0.8, n=9999) +
  geom_point(aes(y = chao1), size = 3, color = 'white') +
  
  labs(y = 'Species richness', x = parse(text = tl)) +
 bktheme + csc +
  yscrich +
  scale_x_continuous(expand = c(0,0), breaks = c(-15,-10,-5,0,5), labels=c(-15,-10,-5,0,5), limits=c(-18.5,8))

panelc <- ggplot(o2015 %>% filter(trait=='logweight'), aes(x=bio6)) + 
  #stat_function(geom='line', fun = fx, args=list(b0 = tempco[1], b1 = tempco[2]), color = 'black', size = 0.8, n=9999) +
  stat_function(geom='line', fun = logitline, args=list(b0 = tempco2[1], b1 = tempco2[2]), color = 'white', size = 0.8, n=9999) +
  geom_point(aes(y = ostat_norm), size = 3, color = 'white') +
 
  labs(y = 'Overlap', x = parse(text = tl)) +
  bktheme + csc +
  scale_x_continuous(expand = c(0,0), breaks = c(-15,-10,-5,0,5), labels=c(-15,-10,-5,0,5), limits=c(-18.5,8))

paneld <- ggplot(o2015 %>% filter(trait=='logweight'), aes(x=ostat_norm)) + 
  stat_function(geom='line', fun = fx2, args=list(b0 = coef2[1], b1 = coef2[2]), color = 'white', size = 0.8, n=9999) +
  geom_point(aes(y = chao1), size = 3, color = 'white') +
  
  labs(x = 'Overlap', y = 'Species richness') +
  bktheme + csc +
  yscrich +
  scale_x_continuous(expand = c(0,0), breaks = c(0, 0.25, 0.5, 0.75), labels=c(0, 0.25, 0.5, 0.75), limits=c(0,0.91))

fp <- 'C:/Users/Q/Dropbox/presentations/sesync2018'
ggsave(file.path(fp, '3panelb.png'), panelb, height=5, width=5, dpi=400)
ggsave(file.path(fp, '3panelc.png'), panelc, height=5, width=5, dpi=400)
ggsave(file.path(fp, '3paneld.png'), paneld, height=5, width=5, dpi=400)


# Map of domains, states, and neon plots ----------------------------------

library(sp)
library(rgdal)
library(ggplot2)
domains <- readOGR(dsn = '/mnt/research/neon/raw_data/spatial_data', layer = 'NEON_Domains')
states <- read.csv('/mnt/research/neon/external_data/raw_external_data/states_albers.csv', stringsAsFactors = FALSE)
rich <- o2015 %>%
  filter(trait == 'logweight') %>%
  dplyr::select(decimalLongitude, decimalLatitude, chao1)
rich <- SpatialPointsDataFrame(coords = as.matrix(cbind(rich$decimalLongitude, rich$decimalLatitude)), 
                         data=data.frame(richness = rich$chao1),
                         proj4string = CRS('+proj=longlat'))


aea_crs <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"
domains <- spTransform(domains, CRSobj = CRS(aea_crs))
rich <- spTransform(rich, CRSobj = CRS(aea_crs))

domains@data <- domains@data %>%
  mutate(id = rownames(domains@data))

domains_fort <- domains %>%
  fortify(region = 'id') %>% 
  left_join(domains@data, by = 'id')

rich_df <- data.frame(long = rich@coords[,1], lat = rich@coords[,2], richness = rich@data$richness)

blktheme <- theme_bw() + 
  theme(axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.title = element_blank(), 
        panel.grid = element_blank(), 
        panel.background = element_rect(color = 'black', fill = 'black'), 
        panel.border = element_blank(), 
        plot.background = element_rect(fill = 'black'), 
        legend.position = c(0.17,0.1), 
        legend.direction = 'horizontal',
        legend.text = element_text(color = 'white'),
        legend.title = element_text(color = 'white'),
        legend.background = element_rect(color = 'black', fill = 'black')
        )

mamcs <- scale_fill_gradientn(name = 'Rodent Species Richness', colours = rev(RColorBrewer::brewer.pal(9, 'RdYlBu')), breaks=c(3,6,9,12))


p <- ggplot(domains_fort) +
  geom_path(data = states %>% filter(!region %in% c('alaska','hawaii')), aes(x = long, y = lat, group = group), color = 'gray50') +
  geom_path(aes(x=long, y=lat, group=group), color = 'white', size = 0.75) +
  geom_point(data = rich_df, aes(x=long, y=lat, fill = richness), color = 'gray40', shape = 21, size = 5, stroke = 2) +
  mamcs +
  coord_equal(xlim = c(-2300000, 2083000), ylim = c(320000, 3100000)) +
  blktheme

ggsave('C:/Users/Q/Dropbox/presentations/sesync2018/neonrichnessmap.png', p, height=6, width=9, dpi=400)
