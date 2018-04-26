# New mammal species richness map.
# 27 Oct 2016
# Modified 02 Dec 2016
# Rerun on 12 Jan 2017

load(file.path(data_path, 'final_data/allorganismal_latest.r'))
load(file.path(data_path, 'final_data/allsiteplot_latest.r'))
richnessdf <- read.csv('richness.csv',stringsAsFactors = FALSE)
richnessdf <- richnessdat

library(ggplot2)

datamap <- ggplot(subset(neonsitedata, decimalLatitude>20 & decimalLatitude<50), aes(x=decimalLongitude, y=decimalLatitude)) +
  borders('state', fill='beige') +
  coord_map() + 
  scale_x_continuous(limits = c(-125,-66.2), expand=c(0,0), breaks = c(-120, -105, -90, -75), labels = c('120° W', '105° W', '90° W', '75° W')) +
  scale_y_continuous(limits=c(24.9,50), expand=c(0,0), breaks = c(25, 35, 45), labels = c('25° N', '35° N', '45° N')) +
  theme(panel.border = element_rect(color='black', fill=NA), panel.background=element_blank(), panel.grid=element_blank(), axis.title = element_blank())



mamcs <- scale_fill_gradientn(name = 'Rodent Species Richness', colours = rev(RColorBrewer::brewer.pal(9, 'RdYlBu')), breaks=c(5,10,15))

neonsitedata <- merge(neonsitedata,richnessdf)

us_map_mr <- datamap + geom_point(size=5, shape=21, aes(fill = chao1), data=subset(neonsitedata, decimalLatitude>20 & decimalLatitude<50 & !is.na(chao1) & !siteID %in% c('DSNY','DELA'))) +
  mamcs +
  theme(legend.position=c(0,0), legend.justification=c(0,0), legend.direction='horizontal')

ggsave('C:/Users/Q/Google Drive/NEON_EAGER/Figures/msfigs2017jan/mapbysite_mammalrichness_newcolors.png', us_map_mr, height=6, width=9, dpi=400)

usmap(filepath = 'C:/Users/Q/Google Drive/NEON_EAGER/Figures/msfigs2017jan/mapbysite_mammalrichness_newcolors.png', fill_var = 'chao1', fill_name = 'Rodent Species Richness', color_all = rev(RColorBrewer::brewer.pal(9, 'RdYlBu')), puertorico = FALSE, alaska = FALSE, inset_titles = FALSE, parse_name = FALSE, legend_bar_extend = 0.5, font = 'Helvetica', mapdat = subset(neonsitedata, !is.na(chao1) & !siteID %in% c('DSNY','DELA')))
