fx <- function(x, b0, b1) 1/(1 + exp(-(b0 + b1 * x)))
fx_loglog <- function(x, b0, b1) exp(-exp(-(b0 + b1 * x)))
tempco <- summary(reglocalbio)$coeff$mean[,1]
chaoco <- summary(reglocalchao)$coeff$mean[,1]
csc <- scale_color_manual(values = c('gray75','black'))

porawtemp <- ggplot(o2015goodsites %>% filter(trait=='logweight'), aes(x=bio1)) + 
  stat_function(geom='line', fun = fx, args=list(b0 = tempco[1], b1 = tempco[2]), color = 'black', size = 0.8, n=9999) +
  geom_point(aes(y = ostat_norm, color=local_significant), size = 3) +
  geom_text(data=data.frame(bio1=20 + 0.5, ostat_norm=0.7, letter='a'), aes(y=ostat_norm,label=letter), size=10) +
  labs(y = 'Overlap', x = parse(text = bioclimnames[1])) +
  theme_john + theme(legend.position = 'none') + csc +
  scale_x_continuous(expand = c(0,0), breaks = c(0,10,20), labels=c(0,10,20), limits=c(-0.5,21.5))

porawchao <- ggplot(o2015goodsites %>% filter(trait=='logweight'), aes(x=chao1)) + 
  stat_function(geom='line', fun = fx, args=list(b0 = chaoco[1], b1 = chaoco[2]), color = 'black', size = 0.8, n=9999) +
  geom_point(aes(y = ostat_norm, color=local_significant), size = 3) +
  geom_text(data=data.frame(chao1=15 + 0.5*(12/22), ostat_norm=0.7, letter='b'), aes(y=ostat_norm,label=letter), size=10) +
  labs(y = 'Overlap', x = 'Species Richness (Chao1)') +
  theme_john + theme(legend.position = 'none') + csc +
  scale_x_continuous(expand = c(0,0), breaks = c(5,10,15), labels=c(5,10,15), limits=c(4,16))

porawtempmanualtransform <- ggplot(o2015goodsites %>% filter(trait=='logweight'), aes(x=bio1)) + 
  geom_abline(intercept = tempco[1], slope = tempco[2], size = 0.8, color = 'black') +
  geom_point(aes(y = qlogis(ostat_norm), color=local_significant), size = 3) +
  #geom_abline(intercept = tempco[1], slope = tempco[2], size = 1.5, color = 'dodgerblue') +
  geom_text(data=data.frame(bio1=20 + 0.5, ostat_norm=0.65, letter='a'), aes(y=qlogis(ostat_norm),label=letter), size=10) +
  scale_y_continuous(breaks = qlogis(c(0.5, 0.1, 0.01)), labels = c(0.5, 0.1, 0.01)) +
  labs(y = 'Overlap', x = parse(text = bioclimnames[1])) +
  theme_john + theme(legend.position = 'none') + csc +
  scale_x_continuous(expand = c(0,0), breaks = c(0,10,20), labels=c(0,10,20), limits=c(-0.5,21.5))

porawchaomanualtransform <- ggplot(o2015goodsites %>% filter(trait=='logweight'), aes(x=chao1)) + 
  geom_abline(intercept = chaoco[1], slope = chaoco[2], size = 0.8, color = 'black') +
  geom_point(aes(y = qlogis(ostat_norm), color=local_significant), size = 3) +
  #geom_abline(intercept = chaoco[1], slope = chaoco[2], size = 1.5, color = 'dodgerblue') +
  geom_text(data=data.frame(chao1=15 + 0.5*(12/22), ostat_norm=0.65, letter='b'), aes(y=qlogis(ostat_norm),label=letter), size=10) +
  scale_y_continuous(breaks = qlogis(c(0.5, 0.1, 0.01)), labels = c(0.5, 0.1, 0.01)) +
  labs(y = 'Overlap', x = 'Species Richness (Chao1)') +
  theme_john + theme(legend.position = 'none') + csc +
  scale_x_continuous(expand = c(0,0), breaks = c(5,10,15), labels=c(5,10,15), limits=c(4,16))


library(gridExtra)
png('figs/msfigs/fig4new.png', height=7, width=4.5, res=400, units='in')
grid.arrange(porawtempmanualtransform, porawchaomanualtransform, nrow=2)
dev.off()

png('figs/msfigs/fig4new_untransformed.png', height=7, width=4.5, res=400, units='in')
grid.arrange(porawtemp, porawchao, nrow=2)
dev.off()
