data {
	int<lower=1> N; // Number of observations
	int n[N];  	 // Sample size
	int y[N];    // outcomes
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
		
	for (i in 1:N){
		lambda1[i] = mu[i]+(1-p)*w*fmin(mu[i]/p,(1-mu[i])/(1-p));
		lambda2[i] = mu[i]-p*w*fmin(mu[i]/p,(1-mu[i])/(1-p));
				
		//con probabilità p
		a1[i] = lambda1[i]*phi;
		b1[i] = (1-lambda1[i])*phi;
		
		//con probabilità 1-p
		a2[i] = lambda2[i]*phi;
		b2[i] = (1-lambda2[i])*phi;
	}
}

model {
//priors
	
	p ~ beta(1,1);
	w ~ beta(1,1);
	
	theta ~ beta(hyper_theta_a, hyper_theta_b);

	for (l in 1:K) {
		if(link_prior_beta == 1)
			beta[l] ~ normal(0, hyperprior_beta);
		else if(link_prior_beta == 2)
			beta[l] ~ cauchy(0, hyperprior_beta);
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
