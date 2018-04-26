# Subset of formatted density shade plots, with colors and with only 4 or 5 plots.
# Edited 19 July to have the minimum rather than mean temperature displayed

sites2use <- c('STEI','BART','KONZ','JORN')
l1 <- labeller(siteID = c(STEI='Steigerwaldt', BART='Bartlett', KONZ='Konza', JORN='Jornada'))
sites_temporder <- neonsitedata %>% arrange(bio6) %>% dplyr::select(siteID, bio6)

pdensshade4 <- ggplot(filter(mam_capture_sitemerge, year == 2015, siteID %in% sites2use) %>% mutate(siteID = factor(siteID, levels=sites_temporder$siteID))) +
  stat_density(adjust = 2, size = 1, aes(x = log10(weight), group = taxonID, fill=taxonID), alpha = 0.25, geom='polygon', position = 'identity') + facet_wrap(~ siteID, ncol = 1, labeller = l1) +
  scale_x_continuous(name = expression(paste('log'[10],' body mass')), breaks = c(1, 2), labels = c(10, 100), limits = c(0.5,2.8)) +
  scale_y_continuous(name = 'probability density', expand = c(0,0), limits=c(0,9)) +
  geom_text(aes(label = paste('NO =', round(ostat_norm,3)), x = 2.5, y = 8), color = 'black', data = o2015 %>% filter(siteID %in% sites2use, trait %in% 'logweight') %>% mutate(siteID = factor(siteID, levels=sites_temporder$siteID)), family = 'Helvetica') +
  geom_text(aes(label = paste0(round(bio1, 2), '°C'), x = 1, y = 8), color = 'black', data = neonsitedata %>% filter(siteID %in% sites2use) %>% mutate(siteID = factor(siteID, levels=sites_temporder$siteID)), family = 'Helvetica') +
  theme_john + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position = 'none', strip.text = element_text(family='Helvetica'))

ggsave('figs/msfigs/4paneldensplot.png', pdensshade4, height = 9, width = 4, dpi = 400)

colorlist <- c('darkorange2', 'gold2', 'black', 'royalblue3','purple3', 'forestgreen', 'red3')
set.seed(27701)
colorvalues <- sample(colorRampPalette(colorlist)(24))

# Edited version with better facet display and location of text on the figure.
pdensshade4clean <- ggplot(filter(mam_capture_sitemerge, year == 2015, siteID %in% sites2use) %>% mutate(siteID = factor(siteID, levels=sites_temporder$siteID))) +
  stat_density(adjust = 2, size = 1, aes(x = log10(weight), group = taxonID, fill=taxonID), alpha = 0.67, geom='polygon', position = 'identity') + facet_wrap(~ siteID, ncol = 1, labeller = l1) +
  scale_fill_manual(values = colorvalues) +
  scale_x_continuous(name = 'Body Mass (g)', breaks = c(1, 2, 3), labels = c(10, 100, 1000), limits = c(0.5,3)) +
  scale_y_continuous(name = 'Probability Density', expand = c(0,0), limits=c(0,9)) +
  geom_text(aes(label = paste('Overlap =', round(ostat_norm,3)), x = 2.5, y = 8.5), color = 'black', data = o2015 %>% filter(siteID %in% sites2use, trait %in% 'logweight') %>% mutate(siteID = factor(siteID, levels=sites_temporder$siteID)), family = 'Helvetica') +
  geom_text(aes(label = paste0('MTCM = ', round(bio6, 1), '°C'), x = 0.5, y = 8.5), color = 'black', data = neonsitedata %>% filter(siteID %in% sites2use) %>% mutate(siteID = factor(siteID, levels=sites_temporder$siteID)), family = 'Helvetica', hjust = 0) +
  theme_john + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position = 'none', strip.text = element_text(family='Helvetica'), strip.background = element_blank())

ggsave('C:/Users/Q/google_drive/NEON_EAGER/Figures/msfigs2017jan/fig3.png', pdensshade4clean, height = 9, width = 4, dpi = 400)


