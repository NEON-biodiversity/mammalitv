# Mediation model
# Source: http://lavaan.ugent.be/tutorial/mediation.html

# exogenous variable (X) is temperature, z-transformed 
# mediator (M) is overlap, logit-transformed
# endogenous variable (Y) is Chao1 richness, z-transformed

mediation_dat <- structure(list(
  X = c(-0.983157854415803, 0.269153855886513, -0.409053263137938, 
        0.282261024562685, -0.549443257222837, 1.61938951573856, 0.933852990010343, 
        0.312012999845949, -0.106897767655547, -2.02603411541746, 0.914298970644232, 
        -0.350863114964422, 0.713144989073109, 1.92892692912945, 0.201262998747074, 
        0.536662194128489, -1.20409661948428, -0.204045158505972, 1.2330578474005, 
        -1.13640708571057, 0.367790398205024, -1.22757605208951, -1.1142404247676
  ), 
  M = c(0.280087359701898, 0.610235589466509, -1.87844875098735, 
        0.638307849044372, 2.13873263027258, -2.95234923234684, -4.01129789291699, 
        -1.24819684323377, -3.13386815593682, 1.30477938055795, -1.71013038540848, 
        -3.8711889786639, -1.73881615699586, -1.33674817486503, -0.874278590275116, 
        -0.168473581540691, 0.778582768874621, -0.886064855035078, -0.193311827816248, 
        1.7700962609674, -1.9238891609905, 0.336906497623836, 0.00114853232086261
  ), 
  Y = c(-2.14807586435342, -2.14807586435342, 0.55397389144382, 
        -2.14807586435342, -1.03546714137809, -1.19441124466028, 0.712917994726011, 
        -0.558634831531514, 1.03080620129039, -1.4063367157032, 0.0771415815972485, 
        -0.240746624967133, -1.83018765778904, -1.19441124466028, -0.876523038095895, 
        -2.4659640709178, -1.51229945122466, -0.558634831531514, -1.83018765778904, 
        -1.83018765778904, -1.83018765778904, -1.19441124466028, -2.14807586435342
  )), 
  class = "data.frame", row.names = c(NA, -23L), .Names = c("X", "M", "Y"))

full_mediation_model <- 
     '# direct effect
      Y ~ c*X
      # mediator
      M ~ a*X
      Y ~ b*M
      # indirect effect (a*b)
      ab := a*b
      # total effect
      total := c + (a*b)
      # intercepts
      M ~ 1
      Y ~ 1
'

library(blavaan)

mediation_bayes <- blavaan(model = full_mediation_model, data = mediation_dat,
                           auto.var = TRUE, auto.fix.first = TRUE, auto.cov.lv.x = TRUE,
                           jagcontrol = list(method = 'rjparallel'),
                           n.chains = 3)

summary(mediation_bayes)
parameterEstimates(mediation_bayes)
inspect(mediation_bayes, 'rsquare')