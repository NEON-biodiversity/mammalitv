# Density plots with different sizes and shadings.

sites2use <- c('STEI','BART','KONZ','JORN')
l1 <- labeller(siteID = c(STEI='Steigerwaldt', BART='Bartlett', KONZ='Konza', JORN='Jornada'))
sites_temporder <- neonsitedata %>% arrange(bio1) %>% dplyr::select(siteID, bio1)

colorlist <- c('darkorange2', 'gold2', 'black', 'royalblue3','purple3', 'forestgreen', 'red3')
set.seed(27701)
colorvalues <- sample(colorRampPalette(colorlist)(26))

dat <- filter(mam_capture_sitemerge, year == 2015, siteID %in% sites2use) %>% mutate(siteID = factor(siteID, levels=sites_temporder$siteID))
abunds <- dat %>% group_by(siteID, taxonID) %>% summarize(n = sum(!is.na(weight)))
total_abunds <- with(abunds, tapply(n,siteID,sum))
abunds <- left_join(abunds, data.frame(siteID=names(total_abunds), total_n=total_abunds))

# Use alpha scale to shade.
pdensshade4alpha <- ggplot(left_join(dat, abunds %>% ungroup %>% mutate(siteID = factor(siteID, levels=sites_temporder$siteID)))) +
  stat_density(adjust = 2, size = 1, aes(x = log10(weight), group = taxonID, alpha=n/total_n), fill = 'blue', geom='polygon', position = 'identity') + 
  facet_wrap(~ siteID, ncol = 1, labeller = l1) +
  scale_fill_manual(values = colorvalues) +
  scale_alpha_continuous(range = c(0.1, 0.75), name = 'relative\nabundance') +
  scale_x_continuous(name = 'body mass (g)', breaks = c(1, 2, 3), labels = c(10, 100, 1000), limits = c(0.5,3.1)) +
  scale_y_continuous(name = 'probability density', expand = c(0,0), limits=c(0,9)) +
  geom_text(aes(label = paste('Overlap =', round(ostat_norm,3)), x = 2.5, y = 8.5), color = 'black', data = o2015 %>% filter(siteID %in% sites2use, trait %in% 'logweight') %>% mutate(siteID = factor(siteID, levels=sites_temporder$siteID)), family = 'Helvetica') +
  geom_text(aes(label = paste0('MAT = ', round(bio1, 1), '°C'), x = 0.5, y = 8.5), color = 'black', data = neonsitedata %>% filter(siteID %in% sites2use) %>% mutate(siteID = factor(siteID, levels=sites_temporder$siteID)), family = 'Helvetica', hjust = 0) +
  theme_john + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position = 'none', strip.text = element_text(family='Helvetica'), strip.background = element_blank())

ggsave('C:/Users/Q/google_drive/NEON_EAGER/Figures/tstat_and_ostat/density_abundance_shading.png', pdensshade4alpha + theme(legend.position = 'bottom'), height = 9, width = 4, dpi = 400)

####################################################

# Same shading, but size of curves depends on abundance

# Must get raw kernel density fns out of the data.

dodens <- function(x) {
  if (sum(!is.na(x$weight)) < 2) return(NA)
  else density(log10(x$weight[!is.na(x$weight)]), adjust=2)
}

library(plyr)

densout <- ddply(dat, .(siteID, taxonID), function(z) {
  dz <- dodens(z)
  if (!inherits(dz,'density')) return(data.frame(x=NA,y=NA))
  else return(data.frame(x=dz$x, y=dz$y))
})

total_abunds <- with(abunds, tapply(n,siteID,sum))
abunds <- left_join(abunds, data.frame(siteID=names(total_abunds), total_n=total_abunds))

pdensshade4height <- ggplot(left_join(densout, abunds %>% mutate(siteID = factor(siteID, levels=sites_temporder$siteID)))) +
  geom_polygon(aes(x=x, y=y*n/total_n, group=taxonID, fill=taxonID), alpha = 0.67) +
  facet_wrap(~ siteID, ncol = 1, labeller = l1) +
  scale_fill_manual(values = colorvalues) +
  scale_x_continuous(name = expression(paste('log'[10],' body mass')), breaks = c(1, 2), labels = c(10, 100), limits = c(0.5,2.7)) +
  scale_y_continuous(name = 'probability density * relative abundance', expand = c(0,0)) +
  geom_text(aes(label = paste('Overlap =', round(ostat_norm,3)), x = 2.35, y = 2.5), color = 'black', data = o2015 %>% filter(siteID %in% sites2use, trait %in% 'logweight') %>% mutate(siteID = factor(siteID, levels=sites_temporder$siteID)), family = 'Helvetica') +
  geom_text(aes(label = paste0('MAT = ', round(bio1, 1), '°C'), x = 0.5, y = 2.5), color = 'black', data = neonsitedata %>% filter(siteID %in% sites2use) %>% mutate(siteID = factor(siteID, levels=sites_temporder$siteID)), family = 'Helvetica', hjust = 0) +
  theme_john + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position = 'none', strip.text = element_text(family='Helvetica'), strip.background = element_blank())

ggsave('C:/Users/Q/google_drive/NEON_EAGER/Figures/tstat_and_ostat/density_abundance_curvesize.png', pdensshade4height, height = 9, width = 4, dpi = 400)


###################################################

# Maintain pretty colors, keeping the saturation and value the same but changing hue. That way alpha can reveal abundance.

hue_colors <- hsv(h = seq(0, 1, length.out = 26), s = 0.8, v = 1)

pdensshade4hue <- ggplot(left_join(dat, abunds %>% ungroup %>% mutate(siteID = factor(siteID, levels=sites_temporder$siteID)))) +
  stat_density(adjust = 2, size = 1, aes(x = log10(weight), group = taxonID, alpha=n/total_n, fill=taxonID), geom='polygon', position = 'identity') + 
  facet_wrap(~ siteID, ncol = 1, labeller = l1) +
  scale_fill_manual(values = hue_colors) +
  scale_alpha_continuous(range = c(0.1, 0.75), name = 'relative\nabundance') +
  scale_x_continuous(name = 'body mass (g)', breaks = c(1, 2, 3), labels = c(10, 100, 1000), limits = c(0.5,3.1)) +
  scale_y_continuous(name = 'probability density', expand = c(0,0), limits=c(0,9)) +
  geom_text(aes(label = paste('Overlap =', round(ostat_norm,3)), x = 2.5, y = 8.5), color = 'black', data = o2015 %>% filter(siteID %in% sites2use, trait %in% 'logweight') %>% mutate(siteID = factor(siteID, levels=sites_temporder$siteID)), family = 'Helvetica') +
  geom_text(aes(label = paste0('MAT = ', round(bio1, 1), '°C'), x = 0.5, y = 8.5), color = 'black', data = neonsitedata %>% filter(siteID %in% sites2use) %>% mutate(siteID = factor(siteID, levels=sites_temporder$siteID)), family = 'Helvetica', hjust = 0) +
  theme_john + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position = 'none', strip.text = element_text(family='Helvetica'), strip.background = element_blank()) +
  theme(legend.position = 'bottom') + guides(fill = FALSE)

ggsave('C:/Users/Q/google_drive/NEON_EAGER/Figures/tstat_and_ostat/density_abundance_differenthues.png', pdensshade4hue, height = 9, width = 4, dpi = 400)

###################################################

# Make the color correspond to abundance, using the same color scheme as in Fig.2.

abund_colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(9, 'RdYlBu')))(26)
fsc1 <- scale_fill_gradientn(name = 'relative\nabundance', colors = rev(RColorBrewer::brewer.pal(9, 'RdYlBu')))
fsc2 <- scale_fill_gradientn(name = 'relative\nabundance', colors = rev(RColorBrewer::brewer.pal(11, 'RdYlGn')))
fsc3 <- scale_fill_gradientn(name = 'relative\nabundance', colors = rev(RColorBrewer::brewer.pal(11, 'Spectral')))


pdensshade4abundcolor <- ggplot(left_join(dat, abunds %>% ungroup %>% mutate(siteID = factor(siteID, levels=sites_temporder$siteID)))) +
  stat_density(adjust = 2, size = 1, aes(x = log10(weight), group = taxonID, fill=n/total_n), alpha=0.67, geom='polygon', position = 'identity') + 
  facet_wrap(~ siteID, ncol = 1, labeller = l1) +
  fsc1 +
  scale_alpha_continuous(range = c(0.1, 0.75), name = 'relative\nabundance') +
  scale_x_continuous(name = 'body mass (g)', breaks = c(1, 2, 3), labels = c(10, 100, 1000), limits = c(0.5,3.1)) +
  scale_y_continuous(name = 'probability density', expand = c(0,0), limits=c(0,9)) +
  geom_text(aes(label = paste('Overlap =', round(ostat_norm,3)), x = 2.5, y = 8.5), color = 'black', data = o2015 %>% filter(siteID %in% sites2use, trait %in% 'logweight') %>% mutate(siteID = factor(siteID, levels=sites_temporder$siteID)), family = 'Helvetica') +
  geom_text(aes(label = paste0('MAT = ', round(bio1, 1), '°C'), x = 0.5, y = 8.5), color = 'black', data = neonsitedata %>% filter(siteID %in% sites2use) %>% mutate(siteID = factor(siteID, levels=sites_temporder$siteID)), family = 'Helvetica', hjust = 0) +
  theme_john + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position = 'none', strip.text = element_text(family='Helvetica'), strip.background = element_blank()) +
  theme(legend.position = 'bottom')

ggsave('C:/Users/Q/google_drive/NEON_EAGER/Figures/tstat_and_ostat/density_abundance_differentcolorsametrans.png', pdensshade4abundcolor, height = 9, width = 4, dpi = 400)
ggsave('C:/Users/Q/google_drive/NEON_EAGER/Figures/tstat_and_ostat/density_abundance_differentcolorsametrans_palette2.png', pdensshade4abundcolor+fsc2, height = 9, width = 4, dpi = 400)
ggsave('C:/Users/Q/google_drive/NEON_EAGER/Figures/tstat_and_ostat/density_abundance_differentcolorsametrans_palette3.png', pdensshade4abundcolor+fsc3, height = 9, width = 4, dpi = 400)