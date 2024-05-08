// Stan model for hermit thrush data

data
{
    int<lower=0> n;   // sample size
    int<lower=0> q;   // number of basis vectors
  
    array[n] int<lower=0, upper=1> y;   // binary outcomes
    vector[n] cover;                    // predictor
    vector[n] elevation;                // predictor
    matrix[n, q] B;                     // basis vectors

    matrix[q, q] V;   // V / tau is the prior covariance matrix for gamma
    vector[q] zero;   // mean vector for the prior on gamma
}

parameters
{
    real beta0;        // intercept
    real beta1;        // slope for cover
    real beta2;        // slope for elevation
    vector[q] gamma;   // basis coefficients

    real<lower=0> tau;   // smoothing parameter for spatial prior
}

transformed parameters
{
    real<lower=0> sigma_sq = inv(tau);   // convert precision to variance
}

model
{
    // logistic regression model

    y ~ bernoulli_logit(beta0 + beta1 * cover + beta2 * elevation + B * gamma);
  
    beta0 ~ normal(0, 1000);   // independent, diffuse Gaussian priors for the
    beta1 ~ normal(0, 1000);   // regression coefficients
    beta2 ~ normal(0, 1000);

    // The prior for gamma is multinormal with mean zero and covariance V / tau.

    gamma ~ multi_normal(zero, sigma_sq * V);
    tau ~ gamma(0.5, 0.0005);
}

generated quantities
{
    // Produce sampled probabilities of hermit thrush presence.

    vector[n] p = inv_logit(beta0 + beta1 * cover + beta2 * elevation + B * gamma);
}
