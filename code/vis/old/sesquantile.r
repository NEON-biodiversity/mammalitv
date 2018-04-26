localses<-subset(o2015goodsites, trait=='logweight')$ostat_norm_localnull_ses
regses<-subset(o2015goodsites, trait=='logweight')$ostat_norm_regnull_ses
contises<-subset(o2015contigoodsites, trait=='logweight')$ostat_norm_regnull_ses
mean(localses)
mean(regses)
mean(contises)
fivenum(localses)
quantile(localses, probs=c(0.025,0.5,0.975))
quantile(regses, probs=c(0.025,0.5,0.975))
quantile(contises, probs=c(0.025,0.5,0.975))
