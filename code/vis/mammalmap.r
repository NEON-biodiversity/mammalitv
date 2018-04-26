# The mammal n code was run in parallel on the cluster.
# Load the coordinates and make into a spatial dataframe.

# Edited 26 Apr 2018: Specify path to HPCC at top
data_path <- '/mnt/research/neon'

mammal_n_all <- list()
for (i in 1:25) {
  load(file = paste0(data_path, '/external_data/raw_external_data/iucn/mammaln', i, '.r'))
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
