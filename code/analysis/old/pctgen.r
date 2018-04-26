# Determine percent generalist of each community

mam_guildtables <- mam_capture %>% 
  group_by(siteID) %>%
  do(guildtable = table(.$Pineda_Main_food, exclude = 'Carnivore'))

mam_guildtables <- with(mam_guildtables, data.frame(siteID=siteID, do.call('rbind', guildtable))) %>%
  group_by(siteID) %>%
  summarize(pctgen = Generalist / (Generalist+Granivore+Herbivore+Insectivore))
  