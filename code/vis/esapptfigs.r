#### FIGURES FOR ESA POWERPOINT
# 24 July 2017

# Edited 26 Apr 2018: Specify path to hpcc at top
data_path <- '/mnt/research/neon'

# New mammal richness map. Use same Chao1 for sites and also put the background richness on there.

rm(list=ls())
source('code/vis/newloadplotdat.r')


crswgs84 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
load(file.path(data_path, 'external_data/raw_external_data/other/mammal_n_all.r'))
mammal_n_spatial <- SpatialPixelsDataFrame(points = mammal_n_all[, c('x','y')], data = mammal_n_all[,3:5], proj4string = crswgs84)


richnessdat <- left_join(richnessdat, o2015 %>% filter(trait=='logweight') %>% dplyr::select(siteID, decimalLatitude, decimalLongitude))


neon_mammaln <- extract(x = brick(mammal_n_spatial), y = with(richnessdat, data.frame(x = decimalLongitude, y = decimalLatitude)))
neon_mammaln[which(richnessdat$siteID=='SERC'),] <- neon_mammaln[which(richnessdat$siteID=='SCBI'),] # Correct this because SERC is too close to coastline.

# Calculate mammal richness, and exclude everything except rodents and shrews
r_rodent <- mam_capture %>% group_by(siteID) %>% filter(family %in% c('Cricetidae', 'Dipodidae', 'Geomyidae', 'Heteromyidae', 'Muridae', 'Sciuridae')) %>% summarize(rodentrichness = length(unique(taxonID)))
r_shrew <- mam_capture %>% group_by(siteID) %>% filter(family %in% c('Soricidae')) %>% summarize(shrewrichness = length(unique(taxonID)))
r_other <- mam_capture %>% group_by(siteID) %>% filter(!family %in% c('Soricidae', 'Cricetidae', 'Dipodidae', 'Geomyidae', 'Heteromyidae', 'Muridae', 'Sciuridae')) %>% summarize(otherrichness = length(unique(taxonID)))

r_all <- r_rodent %>% full_join(r_shrew) %>% full_join(r_other)
r_all[is.na(r_all)] <- 0



neonmapdat <- cbind(richnessdat[,c('siteID','decimalLatitude','decimalLongitude')], neon_mammaln) %>% left_join(r_all)


fillscm <- scale_fill_gradientn(name = 'Rodent richness', colours = rev(RColorBrewer::brewer.pal(9, 'RdYlBu')), limits = range(mammal_n_spatial@data$n_rodents, na.rm=T))

mammal_n_all$n_rodents[mammal_n_all$n_rodents == 0] <- NA

mammalmap <- ggplot(mammal_n_all %>% filter(!is.na(n_rodents)), aes(x=x, y=y)) +
  geom_tile(aes(fill = n_rodents)) +
  borders('state', fill='transparent', colour='black') +
  borders('world', regions = c('mexico','canada'), fill='transparent', colour='black') +
  coord_map(xlim = c(-126, -66), ylim=c(25,50)) + fillscm +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  theme_bw() +
  theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks=element_blank(), legend.position = 'bottom')

neonmapdat$transrodentrichness <- with(neonmapdat, rodentrichness * (max(mammal_n_spatial@data$n_rodents, na.rm=T)/max(rodentrichness, na.rm=T)))


mammalmapwithpoints <- mammalmap + geom_point(data = neonmapdat %>% filter(decimalLatitude > 25 & decimalLatitude < 50, !is.na(rodentrichness)), aes(x=decimalLongitude, y=decimalLatitude, fill=transrodentrichness), shape = 21, size = 5, stroke = 2) + theme(legend.title = element_text(size=16), legend.position = c(0.81, 0.91), legend.direction='horizontal')



mrlm <- lm(rodentrichness ~ n_rodents, data=neonmapdat)
mrrsq <- summary(mrlm)$r.sq

insetplot <- ggplot(neonmapdat, aes(x = n_rodents, y = rodentrichness)) + 
  #xlim(10,50) +
  stat_smooth(method = 'lm', se = FALSE) +
  geom_text(data=data.frame(n_rodents=20, rodentrichness=15), aes(label=label), label = deparse(bquote(R^2 == .(round(mrrsq,2)))), parse=T, size=3) +
  labs(x = 'Grid cell richness', y = 'NEON richness') + geom_point() + theme_bw() + theme(panel.grid=element_blank(), axis.text = element_text(size=5), axis.title = element_text(size=8))

library(grid)

png('rodentmap.png', height=6, width=9, units='in', res=400)

grid.newpage()

v1<-viewport(width = 1, height = 1, x = 0.5, y = 0.5) # Plot area for the main map

print(mammalmapwithpoints, vp=v1) # Print the main map
current.vpTree()
downViewport('panel.6-4-6-4')

print(insetplot, vp=dataViewport(xData=c(-126,-66), yData=c(25,50), clip='off', xscale = c(-78,-68), yscale=c(26,36), x=.89, y=.17, height=0.3, width=0.2))

dev.off()
