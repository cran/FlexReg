data {
	int<lower=1> N; // Number of observations
	array[N] int n;  	 // Sample size
	array[N] int y;    // outcomes
	int<lower=1> K; // Number of covariates
	int<lower=1> H; // number of covariates+intercept for phi
	matrix[N,K]    X; // covariate
	matrix[N,H]    Z; // covariate
	
	int<lower=1> link_code_mu;
	int<lower=1> link_prior_beta;
	real<lower=0> hyperprior_beta;  // Prior standard deviation
	
	int<lower=1> link_code_theta;
	int<lower=1> link_prior_psi;
	real<lower=0> hyperprior_psi;  // Prior standard deviation
}

parameters {
	vector[K] beta;
	vector[H] psi;
	
	real<lower=0,upper=1> w;
	real<lower=0,upper=1> p;
}

transformed parameters {
	vector<lower=0,upper=1>[N]  mu;
	vector<lower=0,upper=1>[N]  lambda1;
	vector<lower=0,upper=1>[N]  lambda2;
	vector<lower=0>[N] b1;
	vector<lower=0>[N] b2;
	vector<lower=0>[N] a1; 
	vector<lower=0>[N] a2;
	vector<lower=0,upper=1>[N] theta;
	vector<lower=0>[N]  phi;
	
	
	if(link_code_mu == 1)
		mu = inv_logit(X * beta);
	else if(link_code_mu == 2)
		mu = Phi(X * beta);
	else if(link_code_mu == 3)
		mu = inv_cloglog(X * beta);
	else if(link_code_mu == 4)
		mu = exp(-exp(X * beta));
		

	if(link_code_theta == 1)
		theta = inv_logit(Z * psi);
	else if(link_code_theta == 2)
		theta = Phi(Z * psi);
	else if(link_code_theta == 3)
		theta = inv_cloglog(Z * psi);
	else if(link_code_theta == 4)
		theta = exp(-exp(Z * psi));
		
		
	for (i in 1:N){
		phi[i] = (1-theta[i])/theta[i];
		
		lambda1[i] = mu[i]+(1-p)*w*fmin(mu[i]/p,(1-mu[i])/(1-p));
		lambda2[i] = mu[i]-p*w*fmin(mu[i]/p,(1-mu[i])/(1-p));
				
		//con probabilità p
		a1[i] = lambda1[i]*phi[i];
		b1[i] = (1-lambda1[i])*phi[i];
		
		//con probabilità 1-p
		a2[i] = lambda2[i]*phi[i];
		b2[i] = (1-lambda2[i])*phi[i];
	}
}

model {
//priors
	
	p ~ beta(1,1);
	w ~ beta(1,1);
	
	for (l in 1:K) {
		if(link_prior_beta == 1)
			beta[l] ~ normal(0, hyperprior_beta);
		else if(link_prior_beta == 2)
			beta[l] ~ cauchy(0, hyperprior_beta);
	}
	
	for (l in 1:H) {
		if(link_prior_psi == 1)
			psi[l] ~ normal(0, hyperprior_psi);
		else if(link_prior_psi == 2)
			psi[l] ~ cauchy(0, hyperprior_psi);
	}
	
// likelihood of log(y)
 for(i in 1:N){
	target += log_mix(p, beta_binomial_lpmf(y[i] | n[i], a1[i], b1[i]),
                         beta_binomial_lpmf(y[i] | n[i], a2[i], b2[i])); 
 }

}  

generated quantities{
	vector[N] log_lik;
	for(i in 1:N){
		log_lik[i] = log_mix(p, beta_binomial_lpmf(y[i] | n[i], a1[i], b1[i]),
								beta_binomial_lpmf(y[i] | n[i], a2[i], b2[i]));
	}
}
