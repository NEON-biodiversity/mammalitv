o2015 <- o2015 %>% rename(siteID=site) %>% left_join(neonsitedata) 


# Set plot parameters
theme_john <- theme_bw() + theme(panel.grid = element_blank(), 
                                 axis.text = element_text(size = 12), 
                                 axis.title = element_text(size = 18),
                                 text = element_text(family = 'Helvetica'))

pdensshade4clean <- ggplot(mam_capture_sitemerge %>% filter(year == 2015, siteID %in% sites2use) %>% mutate(siteID = factor(siteID, levels=sites_temporder$siteID))) +
  stat_density(adjust = 2, size = 1, aes(x = log10(weight), group = taxonID, fill=taxonID), alpha = 0.67, geom='polygon', position = 'identity') + facet_wrap(~ siteID, ncol = 1, labeller = l1) +
  scale_fill_manual(values = colorvalues) +
  scale_x_continuous(name = 'Body Mass (g)', breaks = c(1, 2, 3), labels = c(10, 100, 1000), limits = c(0.5,3)) +
  scale_y_continuous(name = 'Probability Density', expand = c(0,0), limits=c(0,9)) +
  geom_text(aes(label = paste('Overlap =', round(ostat_norm,3)), x = 2.5, y = 8.5), color = 'black', data = o2015 %>% filter(siteID %in% sites2use, trait %in% 'logweight') %>% mutate(siteID = factor(siteID, levels=sites_temporder$siteID)), family = 'Helvetica') +
  geom_text(aes(label = paste0('MAT = ', round(bio1, 1), '°C'), x = 0.5, y = 8.5), color = 'black', data = neonsitedata %>% filter(siteID %in% sites2use) %>% mutate(siteID = factor(siteID, levels=sites_temporder$siteID)), family = 'Helvetica', hjust = 0) +
  theme_john + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position = 'none', strip.text = element_text(family='Helvetica'), strip.background = element_blank())

ggsave('C:/Users/Q/google_drive/NEON_EAGER/Figures/fig3edited.pdf', pdensshade4clean, height = 9, width = 4)

###
# Added 12 March:
# Stealth white on black version of plot
sites2use <- c('STEI','BART','KONZ','JORN')
l1 <- labeller(siteID = c(STEI='Steigerwaldt', BART='Bartlett', KONZ='Konza', JORN='Jornada'))
sites_temporder <- neonsitedata %>% arrange(bio6) %>% dplyr::select(siteID, bio6)

colorlist <- c('darkorange2', 'gold2', 'black', 'royalblue3','purple3', 'forestgreen', 'red3')
set.seed(27701)
colorvalues <- sample(colorRampPalette(colorlist)(24))

theme_john <- theme_bw() + theme(panel.grid = element_blank(), 
                                 axis.text = element_text(size = 12), 
                                 axis.title = element_text(size = 18),
                                 text = element_text(family = 'Helvetica'))

bktheme <- theme(text = element_text(color = 'white'),
      panel.grid = element_blank(),
      panel.border = element_rect(color = 'gray50', fill = 'transparent'),
      panel.background = element_rect(fill = 'black'),
      plot.background = element_rect(fill = 'black', color = 'black'),
      axis.text.x = element_text(size = 12, color = 'white'), 
      axis.title = element_text(size = 18, color = 'white'), 
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      legend.position = 'none',
      strip.background = element_blank(),
      strip.text = element_text(family = 'Helvetica', color = 'white'))

pdensshade4clean <- ggplot(mam_capture_sitemerge %>% filter(year == 2015, siteID %in% sites2use) %>% mutate(siteID = factor(siteID, levels=sites_temporder$siteID))) +
  stat_density(adjust = 2, size = 1, aes(x = log10(weight), group = taxonID, fill=taxonID), alpha = 0.67, geom='polygon', position = 'identity') + facet_wrap(~ siteID, ncol = 1, labeller = l1) +
  scale_fill_manual(values = colorvalues) +
  scale_x_continuous(name = 'Body Mass (g)', breaks = c(1, 2, 3), labels = c(10, 100, 1000), limits = c(0.5,3)) +
  scale_y_continuous(name = 'Probability Density', expand = c(0,0), limits=c(0,9)) +
  geom_text(aes(label = paste('Overlap =', round(ostat_norm,3)), x = 2.5, y = 8.5), color = 'white', data = o2015 %>% filter(siteID %in% sites2use, trait %in% 'logweight') %>% mutate(siteID = factor(siteID, levels=sites_temporder$siteID)), family = 'Helvetica') +
  geom_text(aes(label = paste0('MTCM = ', round(bio6, 1), '°C'), x = 0.5, y = 8.5), color = 'white', data = neonsitedata %>% filter(siteID %in% sites2use) %>% mutate(siteID = factor(siteID, levels=sites_temporder$siteID)), family = 'Helvetica', hjust = 0) +
  bktheme

ggsave('neondensplot.png', pdensshade4clean, height = 9, width = 4, dpi = 400)
