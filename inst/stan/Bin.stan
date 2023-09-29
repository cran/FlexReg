data {
	int<lower=1> N; // Number of observations
	array[N] int n;  	 // Sample size
	array[N] int y;    // outcomes
	int<lower=1> K; // Number of covariates
	matrix[N,K] X;  // Design matrix
	int<lower=1> link_code_mu;
	int<lower=1> link_prior_beta;
	real hyperprior_beta;  // Prior standard deviation

}



parameters {
	vector[K] beta;
}


transformed parameters {
	vector<lower=0,upper=1>[N]  mu;
	
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
	for (l in 1:K) {
		if(link_prior_beta == 1)
			beta[l] ~ normal(0, hyperprior_beta);
		else if(link_prior_beta == 2)
			beta[l] ~ cauchy(0, hyperprior_beta);
	}
	
// likelihood of log(y)
 for(i in 1:N){
	target += y[i] * log(mu[i]) + (n[i] - y[i]) * log(1-mu[i]);
 }

}
  
generated quantities{
	vector[N] log_lik;
	for(i in 1:N){
		log_lik[i] = lgamma(n[i]+1) - lgamma(y[i]+1) - lgamma(n[i]-y[i]+1) + y[i] * log(mu[i]) + (n[i] - y[i]) * log(1-mu[i]);
	}
}
