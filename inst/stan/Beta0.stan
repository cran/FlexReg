data {
	int<lower=1> N; // number of observations
	vector[N]    y; // outcomes
	int<lower=1> K; // number of covariates+intercept for mu
	int<lower=1> K0; // number of covariates+intercept for q1
	matrix[N,K]    X; // covariate
	matrix[N,K0]    X0; // covariate
	int<lower=1> link_code_mu;
	int<lower=1> link_prior_beta;
	int<lower=1> link_prior_omega0;
	real hyperprior_beta;  // Prior standard deviation
	real hyperprior_omega0;  // Prior standard deviation
	int<lower=1> prior_code_phi;
	real<lower=0> hyper_phi;
}

parameters {
	vector[K] beta;
	vector[K0] omega0;
	real<lower=0> phi;
}

transformed parameters {
	vector<lower=0,upper=1>[N]  mu;
	vector<lower=0>[N] a;
	vector<lower=0>[N] b;
	vector<lower=0, upper=1>[N] q0;

	q0 = inv_logit(X0 *omega0);
	
	if(link_code_mu == 1)
		mu = inv_logit(X * beta);
	else if(link_code_mu == 2)
		mu = Phi(X * beta);
	else if(link_code_mu == 3)
		mu = inv_cloglog(X * beta);
	else if(link_code_mu == 4)
		mu = exp(-exp(X * beta));

	for (i in 1:N){
		b[i] = (1-mu[i])*phi;
		a[i] = mu[i]*phi;
	}
}

model {
	//priors
	if(prior_code_phi == 1)
	phi ~ gamma(hyper_phi,hyper_phi);
	else if(prior_code_phi == 2)
	phi ~ uniform(0.0,hyper_phi);

	for (l in 1:K) {
	if(link_prior_beta == 1)
		beta[l] ~ normal(0, hyperprior_beta);
	else if(link_prior_beta == 2)
		beta[l] ~ cauchy(0, hyperprior_beta);
	}

	for (l in 1:K0) {
	if(link_prior_omega0 == 1)
		omega0[l] ~ normal(0, hyperprior_omega0);
	else if(link_prior_omega0 == 2)
		omega0[l] ~ cauchy(0, hyperprior_omega0);
	}
	// likelihood of log(y)

	for(i in 1:N){
	if (y[i]==0)
 	target += log(q0[i]);
 	else
	target += log1m(q0[i])+beta_lpdf(y[i] | a[i],b[i]);
}
}


generated quantities{
	vector[N] log_lik;
	for(i in 1:N){
if (y[i]==0)
 	log_lik[i]= log(q0[i]);
else
		log_lik[i] = log1m(q0[i])+beta_lpdf(y[i] | a[i],b[i]);
	}
}   
