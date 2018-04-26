# Results of additional NEON rodent analysis for revised MS

QDR, created 27 March 2017

## 1. Sensitivity of metric to kernel bandwidth selection method

The reviewers wanted to know whether the overlap metric is sensitive to kernel bandwidth. Clearly, the overlap metric will be sensitive to kernel bandwidth *per se*, since any two distributions can be forced to have any overlap value between 0 and 1 if you choose the correct (unrealistic) bandwidth. So I decided looking sensitivity at bandwidth *per se* was trivial, but on the other hand, it is important to look at whether our metric is sensitive to the method used to select bandwidth. In the original MS, we state that we use the default selection method, `nrd0`, in the R function `density()`. However there are five different methods. I compared the overlap metrics that result from using each of the five methods to select kernel bandwidth for the density functions. Using pairwise correlations between the methods, I showed, at least for the NEON rodent dataset, that the metric is not sensitive to the bandwidth selection method. All correlations were at least 0.96. We can state this in the methods without devoting much space to it.

![kernel method comparison](file:///C:\\Users\\Q\\Google Drive\\NEON_EAGER\\Figures\\tstat_and_ostat\\kernelbandwidthsensitivity.png)  
**Figure 1.** Paired scatterplots of NEON rodent overlap metrics calculated with different bandwidth selection methods below diagonal, and pairwise correlations above diagonal.

## 2. Comparison of pairwise mean difference to pairwise overlap

We stated in the MS that ITV is important, but don't show it. This is because, as we know, the mammals don't have a lot of ITV. The relationships between the drivers and pairwise mean difference are similar to the relationships between the drivers and pairwise overlap. (If ITV is low enough, overlap basically reduces to pairwise mean difference.) We could, however, argue that pairwise overlap is better than difference: if two species are far enough apart, their overlap goes to 0, but their pairwise mean difference keeps on increasing. But this might not be a realistic depiction of their functional difference. If their niches have ceased to overlap, it doesn't matter how far apart the two species are. That could actually be an interesting point.

## 3. New SEM Story

Looking at a simple bivariate scatterplot, we see no relationship between temperature and species richness in this set of mammal communities. We fit an SEM to look at whether there are competing direct and indirect effects. The SEM revealed that:

- Direct pathway: There is no direct effect of temperature on richness; if anything it is negative. 
- Indirect pathway: There is a significantly negative effect of temperature on overlap, and a significantly negative effect of overlap on richness. Therefore, the net indirect effect of temperature on richness, mediated by overlap, is positive (product of two negatives is a positive).
- Net effect: Though the indirect pathway is positive, there is no net effect of temperature on richness (i.e., the model spits back out that there is no net effect if you sum up direct path + indirect path, which we can intuit from seeing the cloud of points on the scatterplot). This might be due to biogeographic quirks or other abiotic drivers. The SEM did allow us to tease apart the positive effect that temperature does have on richness, which in this case is masked by other effects. If we would simply have looked at the bivariate relationship between temperature and richness, we might have concluded that there is no effect.

How we got the net effect: `net = direct + indirect`. Indirect effect is the product of the two legs of the indirect pathway, so `net = direct + indirect1*indirect2`. They basically cancel each other out so the net effect is almost zero; the credible interval extends a long way around zero on either side.

### Formal write-up of SEM methods

We fit a structural equation model in a Bayesian framework to estimate the strengths of direct and indirect overlap-mediated effects of temperature on richness. The Bayesian framework enabled us to estimate credible intervals around each of our parameter estimates. We hypothesized that there is a direct causal pathway between temperature and richness, as well as an indirect pathway consisting of a causal effect of temperature on interspecific body size overlap, and a further causal effect of overlap on richness. The model is specified as follows:
`$R = \beta_{1}T + \beta_{2}O$`  
`$O = \beta_{3}T$`  
The size of the direct effect of temperature on richness corresponds to `$\beta_{1}$`, and the size of the indirect overlap-mediated effect corresponds to the product `$\beta_{2}\beta_{3}$`. The net effect of temperature on richness is the sum of the direct and indirect effects `$\beta_{1} + \beta_{2}beta_{3}$`.   

We examined variance inflation factors and bivariate correlations among potential predictor variables to determine which to include in the final model. (precipitation, seasonality of temperature and precipitation, heterogeneity, primary production, leaf area index). Variables related to precipitation and climate seasonality were highly correlated with other potential predictors, and variables related to topographic and environmental heterogeneity explained very little variation in body size overlap and species richness, so these variables were not included in the SEM. In addition, the productivity variable with the highest bivariate correlations with overlap and richness was yearly maximum leaf area index; we retained this variable and discarded other productivity variables.

### Revised write-up of SEM methods including variable selection and the productivity variable

We fit a structural equation model in a Bayesian framework to estimate the strengths of direct and indirect overlap-mediated effects of temperature on richness. The Bayesian framework enabled us to estimate credible intervals around each of our parameter estimates.  We hypothesized that there are direct causal pathways between temperature and richness and between productivity and richness. The hypothesized indirect pathway between temperature and richness consists of a causal effect of temperature on interspecific body size overlap, and a further causal effect of overlap on richness. We also hypothesized a similar indirect overlap-mediated pathway between productivity and richness. Temperature and productivity had a very low bivariate correlation, so we did not include a pathway from temperature to productivity. The model is specified as follows:

`$R = \beta_{1}T + \beta_{2}P + \beta_{3}O$`  
`$O = \beta_{4}T + \beta_{5}P$`

The size of the direct effect of temperature on richness corresponds to `$\beta_{1}$`, and the size of the indirect overlap-mediated effect corresponds to the product `$\beta_{3}\beta_{4}$`. The net effect of temperature on richness is the sum of the direct and indirect effects `$\beta_{1} + \beta_{3}\beta_{4}$`. Similarly, the net effect of productivity on richness is the sum of the direct and indirect overlap-mediated effects: `$\beta_{2} + \beta_{3}\beta_{5}$`. 

We fit the model in the R package *blavaan* (cite), initializing 3 MCMC chains that each took 25000 samples from the posterior distribution. We estimated each of the `$\beta$` parameters, as well as the derived indirect effect sizes and net effect sizes, at each iteration. The first 5000 samples from each chain were discarded as burn-in. We estimated the 95% credible interval around each parameter estimate, as well as an empirical p-value for each parameter, from the posterior samples. 

### Formal write-up of SEM results

