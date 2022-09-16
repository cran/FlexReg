data {
	int<lower=1> N; // number of observations
	vector[N]    y; // outcomes
	int<lower=1> K; // number of covariates+intercept for mu
	int<lower=1> K1; // number of covariates+intercept for q1
	int<lower=1> K0; // number of covariates+intercept for q0
	int<lower=1> H; // number of covariates+intercept for phi
	matrix[N,K]    X; // covariate for mu
	matrix[N,K1]    X1; // covariate for q1
	matrix[N,K0]    X0; // covariate for q0
	matrix[N,H]    Z; // covariate for phi
	int<lower=1> link_code_mu;
	int<lower=1> link_prior_beta;
	int<lower=1> link_prior_omega1;
	int<lower=1> link_prior_omega0;
	real hyperprior_beta;  // Prior standard deviation
	real hyperprior_omega1;  // Prior standard deviation
	real hyperprior_omega0;  // Prior standard deviation
	int<lower=2> link_code_phi;
	int<lower=1> link_prior_psi;
	real hyperprior_psi;  // Prior standard deviation
}

parameters {
	vector[K] beta;
	vector[K1] omega1;
	vector[K0] omega0;
	vector[H] psi;
	real<lower=0,upper=1> w;
	real<lower=0,upper=1> p;
}


transformed parameters {
	vector<lower=0,upper=1>[N]  mu;
	vector<lower=0>[N]  phi;
	vector<lower=0,upper=1>[N]  lambda1;
	vector<lower=0,upper=1>[N]  lambda2;
	vector<lower=0>[N] b1;
	vector<lower=0>[N] b2;
	vector<lower=0>[N] a1;
	vector<lower=0>[N] a2;
	vector<lower=0>[N]  reg1;
	vector<lower=0>[N]  reg0;
	vector<lower=0>[N] wtilde;
	vector<lower=0, upper=1>[N] q1;
	vector<lower=0, upper=1>[N] q0;

	reg1 =exp(X1 * omega1);
	reg0 =exp(X0 * omega0);
	for (i in 1:N){
	q1[i] = reg1[i]/(1+reg1[i]+reg0[i]);
	q0[i] = reg0[i]/(1+reg1[i]+reg0[i]);
	}
	if(link_code_mu == 1)
		mu = inv_logit(X * beta);
	else if(link_code_mu == 2)
		mu = Phi(X * beta);
	else if(link_code_mu == 3)
		mu = inv_cloglog(X * beta);
	else if(link_code_mu == 4)
		mu = exp(-exp(X * beta));

	if(link_code_phi == 2)
		phi = exp(Z * psi);
	else if(link_code_phi == 3)
		phi = square(Z * psi);

	for (i in 1:N){
	  wtilde[i] = w*fmin(mu[i]/p,(1-mu[i])/(1-p));
	  lambda1[i] = mu[i]+(1-p)*wtilde[i];
	  lambda2[i] = mu[i]-p*wtilde[i];
	//con probabilità p
	b2[i] = (1-lambda1[i])*phi[i];
	a2[i] = lambda1[i]*phi[i];
	//con probabilità 1-p
	b1[i] = (1-lambda2[i])*phi[i];
	a1[i] = lambda2[i]*phi[i];
	}
}

model {
	//priors

	for (l in 1:K) {
	if(link_prior_beta == 1)
		beta[l] ~ normal(0, hyperprior_beta);
	else if(link_prior_beta == 2)
		beta[l] ~ cauchy(0, hyperprior_beta);
	}

	for (l in 1:K1) {
	if(link_prior_omega1 == 1)
		omega1[l] ~ normal(0, hyperprior_omega1);
	else if(link_prior_omega1 == 2)
		omega1[l] ~ cauchy(0, hyperprior_omega1);
	}
	for (l in 1:K0) {
	if(link_prior_omega0 == 1)
		omega0[l] ~ normal(0, hyperprior_omega0);
	else if(link_prior_omega0 == 2)
		omega0[l] ~ cauchy(0, hyperprior_omega0);
	}

	for (l in 1:H) {
	if(link_prior_psi == 1)
		psi[l] ~ normal(0, hyperprior_psi);
	else if(link_prior_psi == 2)
		psi[l] ~ cauchy(0, hyperprior_psi);
	}

	// likelihood of log(y)

for(i in 1:N){
if (y[i]==1)
 
target += log(q1[i]);
 
else{
if (y[i]==0)

target += log(q0[i]);
 
else
target += log(1-q0[i]-q1[i])+log_mix(p, beta_lpdf(y[i] | a2[i],b2[i]),
                beta_lpdf(y[i] | a1[i], b1[i]));
}
} 
}

generated quantities{
	vector[N] log_lik;
for(i in 1:N){
if (y[i]==1)
 
log_lik[i]  = log(q1[i]);
 
else{
if (y[i]==0)
 
log_lik[i]  = log(q0[i]);
 
else
log_lik[i] = log(1-q0[i]-q1[i])+log_mix(p, beta_lpdf(y[i] | a2[i],b2[i]),
                       beta_lpdf(y[i] | a1[i], b1[i]));
}
} 
} 
