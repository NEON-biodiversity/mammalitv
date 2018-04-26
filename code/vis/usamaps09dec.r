# USA Maps for AGU talk

# Edited 26 Apr 2018: Specify path to hpcc at top
data_path <- '/mnt/research/neon'

# US Map function
usmap <- function(filepath, fill_var, fill_name, color_low, color_high, color_all = NULL, alaska = TRUE, puertorico = TRUE, inset_titles = TRUE, parse_name = FALSE, legend_bar_extend = NULL, font = 'Arial', mapdat = NULL) {
  
  siteswithdata <- Reduce(union, list(unique(mam_capture$siteID), unique(bet_pinning$siteID), unique(phe_status$siteID)))
  if (is.null(mapdat)) mapdat <- subset(neonsitedata, siteID %in% siteswithdata)
  
  source(file.path(data_path, 'MS1_RodentOverlap/R_data/qutil.r'))
  
  
  library(directlabels)
  library(ggplot2)
  library(extrafont)
  
  scale_range <- range(neonsitedata[neonsitedata$siteID %in% siteswithdata, fill_var], na.rm=TRUE)
  
  if (parse_name) fill_name <- parse(text = fill_name)
  
  fill_scale <- if(is.null(color_all)) scale_fill_gradient(name=fill_name, low=color_low, high=color_high, limits = c(scale_range[1], scale_range[2])) else scale_fill_gradientn(name = fill_name, colours = color_all, limits = c(scale_range[1], scale_range[2]))
  
  us_map <- ggplot(subset(mapdat, decimalLatitude>20 & decimalLatitude<50), aes(x=decimalLongitude, y=decimalLatitude)) +
    borders('state', fill='beige') +
    geom_point(size=5, shape=21, aes_string(fill = fill_var)) + fill_scale +
    geom_dl(aes(label = siteID, y=decimalLatitude+0.5), method=list('top.bumptwice', 
                                                                    dl.move('UKFS', -93.5, 39),
                                                                    dl.move('NIWO', -105.6, 39),
                                                                    dl.move('TREE', -87.5, 45.3),
                                                                    dl.move('STEI', -92, 45.74),
                                                                    dl.move('BLAN', -80.1, 38.5),
                                                                    dl.move('GRSM', -81, 35.2),
                                                                    dl.move('DELA', -90, 32.54),
                                                                    dl.move('LENO', -88.2, 30.4),
                                                                    fontfamily = font,
                                                                    cex = 1.15,
                                                                    fontface = 'bold')) +
    coord_map() + 
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(limits=c(25,50), expand=c(0,0)) +
    theme(panel.background=element_rect(fill='skyblue'), axis.text = element_blank(), axis.title = element_blank(), axis.ticks=element_blank(), legend.position='bottom', legend.title=element_text(size=24), legend.text=element_text(size=12))
  
  if (!is.null(legend_bar_extend)) us_map <- us_map + theme(legend.key.width = grid::unit(legend_bar_extend, 'inches'))
  
  # Inset map Puerto Rico
  pr_map <- ggplot(subset(mapdat, decimalLatitude<20), aes(x=decimalLongitude, y=decimalLatitude)) +
    borders('world', 'puerto', fill='beige') + fill_scale +
    scale_y_continuous(limits = c(17.9, 18.5)) +
    geom_point(size = 5, shape=21, aes_string(fill = fill_var)) + coord_map() + theme_empty() + theme(panel.border = element_rect(fill='transparent')) +
    geom_dl(aes(label = siteID, y=decimalLatitude+0.5), method=list('top.bumpup', fontfamily = font)) 
  
  if (inset_titles) pr_map <- pr_map + ggtitle('Puerto Rico')
  
  # Inset map Alaska
  ak_map <- ggplot(subset(mapdat, decimalLatitude>50), aes(x=decimalLongitude, y=decimalLatitude)) +
    borders('world', 'USA:Alaska', fill='beige') + fill_scale +
    geom_point(size = 5, shape=21, aes_string(fill = fill_var)) + coord_map() + theme_empty() + theme(panel.border = element_rect(fill='transparent')) +
    geom_dl(aes(label = siteID, y=decimalLatitude+2), method=list('top.bumptwice', fontfamily = font, cex=1.15, fontface='bold', dl.move('DEJU', -141, 65.5), dl.move('HEAL', -163, 63))) +
    xlim(-180,-130)
  
  if (inset_titles) ak_map <- ak_map + ggtitle('Alaska')
  
  us_map <- us_map + theme(text = element_text(family = font))
  ak_map <- ak_map + theme(text = element_text(family = font))
  pr_map <- pr_map + theme(text = element_text(family = font))
  
  
  # Print map with insets
  
  library(grid)
  
  png(filepath, height=6, width=9, units='in', res=400)
  
  grid.newpage()
  
  v1<-viewport(width = 1, height = 1, x = 0.5, y = 0.5) # Plot area for the main map
  
  print(us_map, vp=v1) # Print the main map
  current.vpTree()
  downViewport('panel.3-4-3-4')
  
  # Adjust box size so that it is smaller if there isn't an inset title
  if (inset_titles) {
    ak_y <- 0.15
    ak_h <- 0.3
    ak_w <- 0.3
  }
  else {
    ak_y <- 0.14
    ak_h <- 0.28
    ak_w <- 0.28
  }
  
  # Print the map of Alaska
  if (alaska) print(ak_map, vp=dataViewport(xData=c(-125,-65), yData=c(25,50), clip='off', xscale=c(-180, -130), yscale=c(50,80), x=.09, y=ak_y, height=ak_h, width=ak_w))
  
  # Print the map of Puerto Rico
  if (puertorico) print(pr_map, vp=dataViewport(xData=c(-125,-65), yData=c(25,50), clip='off', xscale=c(-68,-65), yscale=c(17.9,18.6), x=.27, y=.15, height=.15, width=.15))
  
  dev.off()
  
  
}


# Abiotic variable maps ---------------------------------------------------



# Load data.
load(file.path(data_path, 'final_data/allorganismal_latest.r')) # Organismal data
load(file.path(data_path, 'final_data/allsiteplot_latest.r')) # Site covariates

fp <- 'C:/Users/Q/google_drive/NEON_EAGER/Figures/maps/09Dec'

library(ggplot2)
library(extrafont)

source('code/bioclimnames.r')
usmap(filepath = file.path(fp, 'mapbysite_temp.png'), fill_var = 'bio1', fill_name = bioclimnames[1], color_low = 'dodgerblue', color_high = 'indianred', puertorico = FALSE, inset_titles = FALSE, parse_name = TRUE, legend_bar_extend = 0.5, font = 'Century Gothic', mapdat = neonsitedata)
usmap(filepath = file.path(fp, 'mapbysite_precip.png'), fill_var = 'bio12', fill_name = bioclimnames[12], color_low = 'tan', color_high = 'skyblue', puertorico = FALSE, inset_titles = FALSE, parse_name = FALSE, legend_bar_extend = 0.5, font = 'Century Gothic', mapdat = neonsitedata)
usmap(filepath = file.path(fp,'mapbysite_tempseasonality.png'), fill_var = 'bio4', fill_name = bioclimnames[4], color_low = 'forestgreen', color_high = 'lightgoldenrod1', puertorico = FALSE, inset_titles = FALSE, parse_name = FALSE, legend_bar_extend = 0.5, font = 'Century Gothic', mapdat = neonsitedata)
usmap(filepath = file.path(fp, 'mapbysite_tempinterannualcv.png'), fill_var = 'cv_bio1', fill_name = 'Interannual temperature variability (CV)', color_low = 'forestgreen', color_high = 'lightgoldenrod1', puertorico = FALSE, inset_titles = FALSE, parse_name = FALSE, legend_bar_extend = 0.5, font = 'Century Gothic', mapdat = neonsitedata)

spatialstat <- read.csv(file.path(data_path, 'external_data/final_external_data/NEON_spatial_stats_highres.csv'))

neonplotdata <- cbind(neonplotdata, spatialstat)

library(dplyr)

spatialmeans <- neonplotdata %>% group_by(siteID) %>% 
  summarize(ruggedness = mean(tri, na.rm=TRUE), roughness = mean(roughness, na.rm=TRUE))

neonsitedata <- left_join(neonsitedata, spatialmeans)

usmap(filepath = file.path(fp, 'mapbysite_ruggedness.png'), fill_var = 'ruggedness', fill_name = 'Terrain ruggedness index', color_all = RColorBrewer::brewer.pal(9, 'OrRd'), puertorico = FALSE, alaska = TRUE, inset_titles = FALSE, parse_name = FALSE, legend_bar_extend = 0.5, font = 'Century Gothic', mapdat=neonsitedata)


# Data coverage map -------------------------------------------------------

# Find extent of data coverage
mammalsites <- unique(mam_capture$siteID)
beetlesites <- unique(bet_pinning$siteID)
phensites <- unique(phe_status$siteID)
siteswithdata <- Reduce(union, list(mammalsites, beetlesites, phensites))

neonsitedata <- transform(neonsitedata,
                          mammal = siteID %in% mammalsites,
                          beetle = siteID %in% beetlesites,
                          phenology = siteID %in% phensites)



datamap <- ggplot(subset(neonsitedata, decimalLatitude>20 & decimalLatitude<50), aes(x=decimalLongitude, y=decimalLatitude)) +
  borders('state', fill='beige') +
  coord_map() + 
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(limits=c(25,50), expand=c(0,0)) +
  theme(panel.background=element_rect(fill='skyblue'), axis.text = element_blank(), axis.title = element_blank(), axis.ticks=element_blank())

fsc <- scale_fill_manual(values = c('gray80', 'green'))
asc <- scale_alpha_manual(values = c(0.2, 1), labels=c('no', 'yes'))

datamap_all <- datamap + geom_point(aes(alpha = mammal & beetle & phenology), size=5) + labs(alpha = 'Coverage: all three datasets') + asc + theme(legend.position = 'bottom')
datamap_mammal <- datamap + geom_point(aes(alpha = mammal), size=5) + labs(alpha = 'Coverage: mammals') + asc + theme(legend.position = 'bottom')
datamap_beetle <- datamap + geom_point(aes(alpha = beetle), size=5) + labs(alpha = 'Coverage: beetles') + asc + theme(legend.position = 'bottom')
datamap_phen <- datamap + geom_point(aes(alpha = phenology), size=5) + labs(alpha = 'Coverage: phenology') + asc + theme(legend.position = 'bottom')


datamap_all <- datamap + geom_point(aes(alpha = mammal & beetle & phenology), size=5) + ggtitle('All three datasets') + asc + theme(legend.position = c(.92,.1), legend.title = element_blank(), legend.text = element_text(size=16, family='Century Gothic'), plot.title=element_text(size=22, family='Century Gothic'))
datamap_mammal <- datamap + geom_point(aes(alpha = mammal), size=5) + ggtitle('Small mammals') + asc + theme(legend.position = 'none', plot.title=element_text(size=22, family='Century Gothic'))
datamap_beetle <- datamap + geom_point(aes(alpha = beetle), size=5) + ggtitle('Ground beetles') + asc + theme(legend.position = 'none', plot.title=element_text(size=22, family='Century Gothic'))
datamap_phen <- datamap + geom_point(aes(alpha = phenology), size=5) + ggtitle('Plant phenology') + asc + theme(legend.position = 'none', plot.title=element_text(size=22, family='Century Gothic'))


library(gridExtra)

png(file.path(fp, 'mapbysite_datacoverage.png'), height=8, width=15, units='in', res=300)

grid.arrange(datamap_mammal, datamap_beetle, datamap_phen, datamap_all, nrow = 2)

dev.off()


# New mammal latitude plot ------------------------------------------------


date1 <- ymd(mam_capture$date)
date1[is.na(date1)] <- mdy(mam_capture$date[is.na(date1)])
mam_capture$date <- date1

# Mammals have abundance so we use the Chao1 estimator.
# Pool all years for talk.

chao1mammal <- ddply(mam_capture, .(siteID), function(x) {
  xmat <- dcast(x, formula = plotID ~ taxonID)[,-1] # Get rid of date col, use number of rows as abundance
  S_obs <- ncol(xmat)
  f1 <- sum(apply(xmat, 2, sum) == 1)
  f2 <- sum(apply(xmat, 2, sum) == 2)
  return(data.frame(chao1 = S_obs + (f1 * (f1 - 1)) / (2 * (f2 + 1))))
})

chao1mammal <- merge(chao1mammal, neonsitedata)

chao1mamlatlm <- lm(chao1 ~ decimalLatitude, data=subset(chao1mammal, siteID != 'HEAL'))

ggplot(subset(chao1mammal, siteID != 'HEAL'), aes(x=decimalLatitude, y=chao1)) +
  geom_line(aes(y = predict(chao1mamlatlm)), size=1.5, color='slateblue') +
  geom_point(size=3) + #geom_dl(aes(label=siteID), method='top.bumptwice') +
  theme_bw() +
  theme(text = element_text(family = 'Century Gothic'), axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), axis.title.x = element_text(size=14), axis.title.y=element_text(size=14)) +
  labs(x = 'Latitude', y = 'Chao1 Richness Estimator') +
  geom_text(data=data.frame(decimalLatitude = 37.5, chao1 = 17, lab = "list(R^2 == 0.17, p == 0.047)"), aes(label=lab), parse=TRUE, family = 'Century Gothic') 
#  ggtitle('Mammal richness versus latitude')
ggsave(file.path(fp, 'latitude_mammalrichness.png'), height=5, width=5, dpi=300)


# New plant latitude plot -------------------------------------------------

chao2 <- ddply(subset(ppc_400m2, year(date) == 2014), .(siteID), function(x) {
  xmat <- dcast(x, formula = date ~ taxonID, value.var='targetTaxaPresent')[,-1] > 0 # Get rid of date col, ignore entries greater than 1
  S_obs <- ncol(xmat)
  m <- nrow(xmat)
  n_sites <- apply(xmat, 2, sum)
  q1 <- sum(n_sites == 1)
  q2 <- sum(n_sites == 2)
  return(data.frame(S_obs = S_obs, chao2 = S_obs + (m - 1) / m * (q1 * (q1 - 1)) / (q2 * (q2 + 1))))
})

chao2 <- merge(chao2, neonsitedata)

chao2latlm <- lm(chao2 ~ decimalLatitude, data=chao2)

source(file.path(data_path, 'MS1_RodentOverlap/R_data/qutil.r'))
ggplot(chao2, aes(x=decimalLatitude, y=chao2)) +
  geom_line(aes(y = predict(chao2latlm)), size=1.5, color='slateblue') +
  geom_point(size=3) + #geom_dl(aes(label=siteID), method='top.bumptwice') +
  theme_bw() +
  theme(text = element_text(family = 'Century Gothic'), axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), axis.title.x = element_text(size=14), axis.title.y=element_text(size=14)) +
  labs(x = 'Latitude', y = 'Chao2 Richness Estimator') +
  geom_text(data=data.frame(decimalLatitude = 39, chao2 = 310, lab = "list(R^2 == 0.39, p == 0.02)"), aes(label=lab), parse=TRUE, family = 'Century Gothic') 
#  ggtitle(expression(paste('Plant richness in 400-m'^2,' plots versus latitude')))
ggsave(file.path(fp,'latitude_plantrichness.png'), height=5, width=5, dpi=300)


# New beetle latitude plot ------------------------------------------------

bet_pinning <- subset(bet_pinning, !is.na(siteID) & siteID!='')

# Correct dates that are not on the same bout in the bet_pinning data frame.

bet_pinning$collectDate[bet_pinning$siteID=='OSBS' & bet_pinning$collectDate=='2014-06-07'] <- '2014-06-09'
bet_pinning$collectDate[bet_pinning$siteID=='DSNY' & bet_pinning$collectDate=='2014-06-06'] <- '2014-06-11'
bet_pinning$collectDate[bet_pinning$siteID=='JERC' & bet_pinning$collectDate=='2014-09-10'] <- '2014-09-11'
bet_pinning$collectDate[bet_pinning$siteID=='ORNL' & bet_pinning$collectDate=='2014-06-17'] <- '2014-06-16'
bet_pinning$collectDate[bet_pinning$siteID=='SCBI' & bet_pinning$collectDate=='2014-06-11'] <- '2014-06-10'
bet_pinning$collectDate[bet_pinning$siteID=='SCBI' & bet_pinning$collectDate=='2014-08-21'] <- '2014-08-19'
bet_pinning$collectDate[bet_pinning$siteID=='UNDE' & bet_pinning$collectDate=='2014-08-12'] <- '2014-08-19'
bet_pinning$collectDate[bet_pinning$siteID=='WOOD' & bet_pinning$collectDate=='2014-08-17'] <- '2014-08-13'
bet_pinning$collectDate[bet_pinning$siteID=='WOOD' & bet_pinning$collectDate=='2014-08-23'] <- '2014-08-27'



chao1beetle <- ddply(bet_pinning, .(siteID), function(x) {
  xmat <- dcast(x, formula = collectDate ~ taxonID)[,-(1:2)] # Get rid of date and var2 col, use number of rows as abundance
  S_obs <- ncol(xmat)
  f1 <- sum(apply(xmat, 2, sum) == 1)
  f2 <- sum(apply(xmat, 2, sum) == 2)
  return(data.frame(chao1 = S_obs + (f1 * (f1 - 1)) / (2 * (f2 + 1))))
})

chao1beetle <- merge(subset(chao1beetle,chao1 > 0), neonsitedata) # get rid of those with zero

chao1betlatlm <- lm(chao1 ~ decimalLatitude, data=chao1beetle)
chao1betlm2 <- lm(chao1 ~ decimalLatitude + I(decimalLatitude^2), data=chao1beetle)

newx <- data.frame(decimalLatitude=seq(27.5,47.5,by=0.1))
ggplot(chao1beetle, aes(x=decimalLatitude, y=chao1)) +
  geom_line(aes(y = predict(chao1betlm2, newdata = newx)), data=newx, size=1.5, color='slateblue') +
  geom_point(size=3) + #geom_dl(aes(label=siteID), method='top.bumptwice') +
  theme_bw() +
  theme(text = element_text(family = 'Century Gothic'), axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), axis.title.x = element_text(size=14), axis.title.y=element_text(size=14)) +
  labs(x = 'Latitude', y = 'Chao1 Richness Estimator') +
  geom_text(data=data.frame(decimalLatitude = 37.5, chao1 = 52, lab = "list(R^2 == 0.60, p == 0.03)"), aes(label=lab), parse=TRUE, family = 'Century Gothic') 
  #qSubtitle('Beetle richness versus latitude', 'not including morphospecies')
ggsave(file.path(fp,'latitude_beetlerichness.png'), height=5, width=5, dpi=300)
# New number of records datatable -----------------------------------------

betrec <- table(bet_pinning$siteID)
mamrec <- table(mam_capture$siteID)
hbprec <- table(hbp_massdata$siteID)
#table(mpr_perrootsample$siteID)
pherec <- table(phe_status$siteID)
ppcrec <- table(ppc_1m2$siteID)

# Create matrix
allsites <- Reduce(union, list(names(betrec), names(mamrec), names(hbprec), names(pherec), names(ppcrec)))
numrecords <- do.call('cbind', list(betrec[allsites[-1]], mamrec[allsites[-1]], hbprec[allsites[-1]], pherec[allsites[-1]], ppcrec[allsites[-1]]))
dimnames(numrecords)[[1]] <- allsites[-1]
dimnames(numrecords)[[2]] <- c('beetle','mammal','clip_harvest','plant_phenology','plant_community')
numrecords[is.na(numrecords)] <- 0

numrecords <- cbind(allsites[-1], numrecords)

# Sort by latitude
northtosouth <- neonsitedata$siteID[order(neonsitedata$decimalLatitude, decreasing = T)]

numrecords <- numrecords[northtosouth[northtosouth %in% numrecords[,1]],]

library(XLConnect)
writeWorksheetToFile(data=numrecords, file=file.path(fp, 'datatable.xlsx'), sheet='Sheet1')
