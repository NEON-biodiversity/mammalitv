// STAN model specification for t-stats multiple regression
// Reduced model for TIPIC by weight

data {
	int N;					// Number of plot-level measurements
	int Nsites;				// Number of distinct site IDs
	int Nparams;			// Number of parameters in model
	
	// Response variable
	real t[N];				// T-statistic (any type)
	
	// Random effect
	int siteID[N];			// Integer representing site ID for each plot-level measurement
	
	// Fixed predictors. These should all be z-transformed.
	real bio1[N];			// Mean annual temperature at each plot
	real bio4[N];			// Within-year temperature variability at each plot
	real cv_bio1[N];		// Among-year temperature variability at each plot
	real bio12[N];			// Mean annual precipitation at each plot
	real bio15[N];			// Within-year precipitation variability at each plot
	real cv_bio12[N];		// Among-year precipitation variability at each plot
	real chao1[N];			// Chao1 species richness estimator for each plot
	real mpd_z[N];			// Phylogenetic MPD index (z-score) for each plot
	real mntd_z[N];			// Phylogenetic MNTD index (z-score) for each plot
	real ruggedness[N];		// Terrain ruggedness index for each plot
	
}

parameters {
	real alpha[Nsites];		// Random intercepts for sites
	real beta[Nparams];		// Regression coefficients (slopes for each parameter, no intercept here)
	real<lower=0> sigma;	// Standard deviation
	
	
}

transformed parameters {
	real mu[N];
	
	for (i in 1:N) {
		// This is the line edited to remove parameters. I left the Nparams the same.
		mu[i] = alpha[siteID[i]] + beta[1] * bio1[i] + beta[2] * bio4[i] + beta[3] * cv_bio1[i] + beta[5] * bio15[i] + beta[6] * cv_bio12[i] + beta[7] * chao1[i] + beta[8] * mpd_z[i] + beta[10] * ruggedness[i];
	}
}

model {
	
	// Likelihood
	t ~ normal(mu, sigma);	// Plain-vanilla regression
	
	// Priors
	// Currently left flat (bounded at zero for the standard deviation).
	
	// Hyperpriors
	// Not needed
	
}

generated quantities {
	// Log-likelihood (used for calculating information criteria).
	real log_lik;
	log_lik = normal_lpdf(t | mu, sigma);
}

