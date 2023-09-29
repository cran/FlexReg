data {
	int<lower=1> N; // Number of observations
	array[N] int n;  	 // Sample size
	array[N] int y;    // outcomes
	int<lower=1> K; // Number of covariates
	matrix[N,K] X;  // Design matrix
	
	int<lower=1> link_code_mu;
	int<lower=1> link_prior_beta;
	real<lower=0> hyperprior_beta;  // Prior standard deviation
	
	real<lower=0> hyper_theta_a;
	real<lower=0> hyper_theta_b;
}



parameters {
	vector[K] beta;
	real<lower=0,upper=1> theta;
}


transformed parameters {
	vector<lower=0,upper=1>[N]  mu;
	real<lower=0> phi; 
	
	phi = (1-theta)/theta;
	
	if(link_code_mu == 1)
		mu = inv_logit(X * beta);
	else if(link_code_mu == 2)
		mu = Phi(X * beta);
	else if(link_code_mu == 3)
		mu = inv_cloglog(X * beta);
	else if(link_code_mu == 4)
		mu = exp(-exp(X * beta));
		
}

model {
//priors
	theta ~ beta(hyper_theta_a, hyper_theta_b);

	for (l in 1:K) {
		if(link_prior_beta == 1)
			beta[l] ~ normal(0, hyperprior_beta);
		else if(link_prior_beta == 2)
			beta[l] ~ cauchy(0, hyperprior_beta);
	}
	
// likelihood of log(y)
 for(i in 1:N){
	target += lgamma(phi) - lgamma(phi * mu[i]) - lgamma(phi * (1-mu[i])) + lgamma(phi * mu[i]+y[i]) + lgamma(phi * (1-mu[i]) + n[i]-y[i]) - lgamma(phi+n[i]);
 }
}

generated quantities{
	vector[N] log_lik;
	for(i in 1:N){
		log_lik[i] = lgamma(n[i]+1) - lgamma(y[i]+1) - lgamma(n[i]-y[i]+1) + lgamma(phi) - lgamma(phi * mu[i]) - lgamma(phi * (1-mu[i])) + lgamma(phi * mu[i]+y[i]) + lgamma(phi * (1-mu[i]) + n[i]-y[i]) - lgamma(phi+n[i]);
	}
}
