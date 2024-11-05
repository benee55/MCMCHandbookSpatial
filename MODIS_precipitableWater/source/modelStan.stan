data {
  int N; //the number of observations
  int K; //the number of columns in the model matrix
  int P; //the number of columns in the basis matrix
  matrix[N,K] X; //the model matrix
  matrix[N,P] M; //the basis matrix
  real<lower=0> y[N]; //the response
}

parameters {
  vector[K] betas; //the regression parameters
  vector[P] deltas; //the reparameterized random effects
  real<lower=0>  tau; //the precision parameter
  real<lower=0> inverse_phi; //the dispersion parameter
}

transformed parameters {
  vector[N] linpred;
  vector[N] mu; //the expected values (linear predictor)
  vector[N] beta; //rate parameter for the gamma distribution
  linpred = X*betas+M*deltas;
  mu = exp(linpred); //using the log link 
  beta = rep_vector(inverse_phi, N)./ mu; // Compound Elementwise Division
}

model {  
  betas[1:] ~ normal(0,100); //prior for the betas following Gelman 2008
  tau ~ gamma(0.5,2000); //prior for the precision parameter for tau2 for new spatial random effects (Hughes and Haram, 2012)
  deltas[1:] ~ normal(0, 1/tau); // MVN prior for the new spatial random effects
  inverse_phi ~ exponential(1); // prior on inverse dispersion parameter
  y ~ gamma(inverse_phi,beta);
}
