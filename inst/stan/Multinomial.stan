data{
	int<lower=1> N; // total number of observations
	//int<lower=2> D; // number of categories
	// D = 3
	int<lower=1> K0; // number of predictor levels
	int<lower=1> K1; // number of predictor levels

	matrix[N,K0] X0; // predictor design matrix
	matrix[N,K1] X1; // predictor design matrix

	//matrix[N,D] Y; // response variable
	int Y[N,3];

	real hyperprior_omega1;
	real hyperprior_omega0;

	int<lower=1> link_prior_omega1;
	int<lower=1> link_prior_omega0;
}

parameters {
	vector[K0] omega0; // coefficients
	vector[K1] omega1; // coefficients
}

transformed parameters{
	//real exptheta = exp(theta);
	vector<lower=0>[N] reg1;
	vector<lower=0>[N] reg0;

	//vector<lower=0, upper=1>[N] q1;
	//vector<lower=0, upper=1>[N] q0;

	matrix[N,3] q;

	reg1 = exp(X1 * omega1);
	reg0 = exp(X0 * omega0);

	for(i in 1:N){
	  q[i,1] = reg0[i]/(1+ reg1[i] + reg0[i]);
	  q[i,2] = reg1[i]/(1+ reg1[i] + reg0[i]);
	  q[i,3] = 1/(1+ reg1[i] + reg0[i]);
	}

}

model {
// prior:
	 if(link_prior_omega1 == 1)
		omega1 ~ normal(0, hyperprior_omega1);
	 else if(link_prior_omega1 == 2)
		omega1 ~ cauchy(0, hyperprior_omega1);


	 if(link_prior_omega0 == 1)
		omega0 ~ normal(0, hyperprior_omega0);
	 else if(link_prior_omega0 == 2)
		omega0 ~ cauchy(0, hyperprior_omega0);




// likelihood
	for (n in 1:N) {
		//transpose(Y[n,]) ~ dirichlet(transpose(mu[n,]) * exptheta);
		Y[n,] ~ multinomial(transpose(q[n,]));
	}
}

generated quantities{
	vector[N] log_lik;
	for(n in 1:N){
		//log_lik[n] = dirichlet_lpdf(transpose(Y[n,]) | transpose(mu[n,]) * exptheta);
		log_lik[n] = multinomial_lpmf(Y[n,] | transpose(q[n,]));
	}
}
