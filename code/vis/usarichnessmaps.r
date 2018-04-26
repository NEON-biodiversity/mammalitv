# Maps with diversity of different taxa

# Edited 26 Apr 2018: Specify path to hpcc at top
data_path <- '/mnt/research/neon'

# Load plant shapefile

library(maptools)
ellis <- readShapeSpatial(fn = file.path(data_path, 'external_data/raw_external_data/other/ellis_2012_shapefile/ellis_2012_l8_dataset_2012_01_17'))

# Contains 16805 hexagons
spplot(ellis, zcol = 'N')

# Mask to the borders of the USA
library(maps)
data(usaMapEnv)

# Get extent of USA
usa_lat <- c(25,50)
usa_long <- c(-125, -67)

library(raster)
usa_extent <- extent(matrix(c(usa_long, usa_lat), nrow=2, byrow=T))

ellis_usa <- crop(x = ellis, y = usa_extent)
spplot(ellis_usa, zcol = 'N', col = NA)

# Create ggplot object
library(ggplot2)
ellis_fort <- fortify(ellis_usa, region = 'L8_ID')

# Merge with additional attributes of ellis
library(dplyr)
ellis_dat <- ellis_usa@data %>% rename(id = L8_ID) %>% mutate(id = as.character(id))
ellis_join <- left_join(ellis_fort, ellis_dat, by = 'id')

# ggplot of ellis, with USA borders on it.
ellis_nrange <- range(ellis_join$N)
fillsc <- scale_fill_gradientn(name = 'Plant richness', colours = rev(RColorBrewer::brewer.pal(9, 'RdYlBu')), limits = ellis_nrange)

plantmap <- ggplot(ellis_join, aes(x=long, y=lat)) +
  geom_polygon(aes(fill = N, group = group)) +
  borders('state', fill='transparent', colour='black') +
  coord_map() + fillsc +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  theme(panel.background=element_rect(fill='gray75'), axis.text = element_blank(), axis.title = element_blank(), axis.ticks=element_blank(), legend.position = 'bottom')

# Get richness for the neon sites, along with their coordinates, and also compare this to the background richness.
load(file.path(data_path, 'final_data/allsiteplot_latest.r'))
load(file.path(data_path, 'final_data/allorganismal_latest.r'))

# Background richness from ellis
neon_ellis <- extract(x = ellis, y = with(neonsitedata, data.frame(x = decimalLongitude, y = decimalLatitude)))

# Plant richness at 400m2 plot scale
# Plants
r400 <- ppc_400m2 %>% group_by(siteID) %>% summarize(richness = length(unique(taxonID)))
neonsitedata <- neonsitedata %>% cbind(neon_ellis) %>% left_join(r400)

plot(richness ~ N, data=neonsitedata)

colorsc <- scale_colour_gradientn(name = 'Plant richness', colours = rev(RColorBrewer::brewer.pal(9, 'RdYlBu')), limits = range(neonsitedata$richness, na.rm=T))

# Translate the local richness to the Ellis scale.
neonsitedata <- mutate(neonsitedata, transrichness = richness * (ellis_nrange[2]/max(richness, na.rm=T)))

plantmapwithpoints <- plantmap + geom_point(data = neonsitedata %>% filter(decimalLatitude > 25 & decimalLatitude < 50, !is.na(richness)), aes(x=decimalLongitude, y=decimalLatitude, fill=transrichness-1), shape = 21, size = 5, stroke = 2) + theme(legend.title = element_text(size=16), legend.position = c(0.81, 0.91), legend.direction='horizontal')

# Inset plot local richness vs grid cell richness

rlm <- lm(richness~N, data=neonsitedata)
rrsq <- summary(rlm)$r.sq

insetplot <- ggplot(neonsitedata, aes(x = N, y = richness)) + 
  stat_smooth(method = 'lm', se = FALSE) + xlim(500,1500) +
  geom_text(data=data.frame(N=750, richness=290), aes(label=label), label = deparse(bquote(R^2 == .(round(rrsq,2)))), parse=T, size=3) +
  labs(x = 'Grid cell richness', y = 'NEON richness') + geom_point() + theme_bw() + theme(panel.grid=element_blank(), axis.text = element_text(size=5), axis.title = element_text(size=8))


library(grid)

png('figs/png/plantmapwithinset.png', height=6, width=9, units='in', res=400)

grid.newpage()

v1<-viewport(width = 1, height = 1, x = 0.5, y = 0.5) # Plot area for the main map

print(plantmapwithpoints, vp=v1) # Print the main map
current.vpTree()
downViewport('panel.3-4-3-4')

print(insetplot, vp=dataViewport(xData=c(-125,-67), yData=c(25,50), clip='off', xscale = c(-78,-68), yscale=c(26,36), x=.89, y=.17, height=0.3, width=0.2))

dev.off()


# Mammal map using iucn ---------------------------------------------------


library(maptools)
library(rgdal)

# Define projection and read the IUCN polygons in
crswgs84 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
x <- readShapePoly(file.path(data_path, 'external_data/raw_external_data/iucn/TERRESTRIAL_MAMMALS.shp'), proj4string = crswgs84)

# Subset the mammal polygons to just rodents and shrews
x <- subset(x, order_ %in% c('RODENTIA', 'EULIPOTYPHLA'))

# Loop through the coordinates of the USA by 1x1 degree boxes and see how many IUCN ranges overlap each box
coords <-  expand.grid(x = (-124.5):(-67.5), y = 25.5:49.5)
xp <- as(x, 'SpatialPolygons')

mammal_n <- rep(0, length(coords))
pb <- txtProgressBar(min=0, max=nrow(coords), style=3)

for (i in 1:nrow(coords)) {
  point_i <- SpatialPoints(cbind(x = coords[i,1], y = coords[i,2]), proj4string = crswgs84)
  n_i <- 0
  for (j in 1:length(xp)) {
    if (!is.na(over(point_i, xp[j]))) n_i <- n_i + 1
  }
  mammal_n[i] <- n_i
  setTxtProgressBar(pb, i)
}

close(pb)

# The mammal n code was run in parallel on the cluster.
# Load the coordinates and make into a spatial dataframe.

mammal_n_all <- list()
for (i in 1:25) {
  load(file = paste0(data_path, 'external_data/raw_external_data/iucn/mammaln', i, '.r'))
  mammal_n_all[[i]] <- mammal_n_dat
}
mammal_n_all <- do.call('rbind', mammal_n_all)
save(mammal_n_all, file = file.path(data_path, 'external_data/raw_external_data/iucn/mammal_n_all.r'))

load(file.path(data_path, 'external_data/raw_external_data/iucn/mammal_n_all.r'))
mammal_n_spatial <- SpatialPixelsDataFrame(points = mammal_n_all[, c('x','y')], data = mammal_n_all[,3:5], proj4string = crswgs84)

neon_mammaln <- extract(x = brick(mammal_n_spatial), y = with(neonsitedata, data.frame(x = decimalLongitude, y = decimalLatitude)))
neon_mammaln[27,] <- neon_mammaln[26,] # Correct this because SERC is too close to coastline.

# Calculate mammal richness, and exclude everything except rodents and shrews
r_rodent <- mam_capture %>% group_by(siteID) %>% filter(family %in% c('Cricetidae', 'Dipodidae', 'Geomyidae', 'Heteromyidae', 'Muridae', 'Sciuridae')) %>% summarize(rodentrichness = length(unique(taxonID)))
r_shrew <- mam_capture %>% group_by(siteID) %>% filter(family %in% c('Soricidae')) %>% summarize(shrewrichness = length(unique(taxonID)))
r_other <- mam_capture %>% group_by(siteID) %>% filter(!family %in% c('Soricidae', 'Cricetidae', 'Dipodidae', 'Geomyidae', 'Heteromyidae', 'Muridae', 'Sciuridae')) %>% summarize(otherrichness = length(unique(taxonID)))

r_all <- r_rodent %>% full_join(r_shrew) %>% full_join(r_other)
r_all[is.na(r_all)] <- 0

#neonsitedata <- neonsitedata %>% mutate(mammal_gridcell = neon_mammaln) %>% left_join(r_mammal %>% rename(mammal_richness=richness))

#plot(mammal_richness ~ mammal_gridcell, data=neonsitedata)

neonmapdat <- cbind(neonsitedata[,c('siteID','decimalLatitude','decimalLongitude')], neon_mammaln) %>% left_join(r_all)

fillscm <- scale_fill_gradientn(name = 'Rodent richness', colours = rev(RColorBrewer::brewer.pal(9, 'RdYlBu')), limits = range(mammal_n_spatial@data$n_rodents, na.rm=T))

mammal_n_all$n_rodents[mammal_n_all$n_rodents == 0] <- NA

mammalmap <- ggplot(mammal_n_all %>% filter(!is.na(n_rodents)), aes(x=x, y=y)) +
  geom_tile(aes(fill = n_rodents)) +
  borders('state', fill='transparent', colour='black') +
  coord_map() + fillscm +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  theme(panel.background=element_rect(fill='gray75'), axis.text = element_blank(), axis.title = element_blank(), axis.ticks=element_blank(), legend.position = 'bottom')

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

png('figs/png/mammalmapwithinset_rodents_30minres.png', height=6, width=9, units='in', res=400)

grid.newpage()

v1<-viewport(width = 1, height = 1, x = 0.5, y = 0.5) # Plot area for the main map

print(mammalmapwithpoints, vp=v1) # Print the main map
current.vpTree()
downViewport('panel.3-4-3-4')

print(insetplot, vp=dataViewport(xData=c(-125,-67), yData=c(25,50), clip='off', xscale = c(-78,-68), yscale=c(26,36), x=.89, y=.17, height=0.3, width=0.2))

dev.off()