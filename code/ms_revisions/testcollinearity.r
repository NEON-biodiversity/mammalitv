# Check collinearity of predictors and do other variable selection (to show we only use temperature to predict richness)
# Update 03 April: add longitude as a predictor (proximity to Sierra Madre Corridor)


source('code/vis/newloadplotdat.r')
data_path <- '/mnt/research/neon'

# add productivity means.

modisdat <- read.csv(file.path(data_path, 'external_data/final_external_data/NEON_modisyearly_withcv.csv'), stringsAsFactors = FALSE)
modis_mean <- modisdat %>% group_by(siteID) %>% summarize_at(-(1:3), mean, na.rm=TRUE)

all_preds <- o2015 %>%
  filter(trait == 'logweight') %>%
  left_join(modis_mean) %>%
  dplyr::select(elevation, bio1, bio4, bio6, bio12, bio15, bio14, cv_bio1, cv_bio12, pc1_productivityheterogeneity, pc2_topographyheterogeneity, LAI, NDVI, NPP, GPP, cv_LAI, cv_NDVI, decimalLongitude)

library(usdm)
vif(all_preds)

# Remove NPP and GPP.
all_preds <- dplyr::select(all_preds, -NPP, -GPP)
vif(all_preds)

all_preds <- dplyr::select(all_preds, -elevation, -bio1)
vif(all_preds)

all_preds <- dplyr::select(all_preds, -cv_LAI)
vif(all_preds)

all_preds <- dplyr::select(all_preds, -bio14)
vif(all_preds)

all_preds <- dplyr::select(all_preds, -bio4)
vif(all_preds)

all_preds <- dplyr::select(all_preds, -bio12, -NDVI)
vif(all_preds)


library(GGally)

all_variables <- o2015 %>%
  filter(trait == 'logweight') %>%
  left_join(modis_mean) %>%
  dplyr::select(bio6, bio15, cv_bio1, cv_bio12, pc1_productivityheterogeneity, pc2_topographyheterogeneity, decimalLongitude, NPP, cv_NDVI, chao1, ostat_norm)

ggpairs(all_variables)

cor(all_variables)

# It's pretty subjective at this point what variable(s) to use. We have temperature and productivity which correlates with precip. So it's probably OK to include temperature and productivity.

# Variables that survived the multicollinearity test and have strong bivariate correlations are:
# mean annual temperature and yearly maximal leaf area index. The second is an index of productivity. It essentially changes along the longitudinal gradient.



