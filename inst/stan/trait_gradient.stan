// Bivariate phylogenetic model of the trait gradient (comdynamite::trait_gradient_bayes).
// Response 1 is the site mean of a population (without the focal species), response 2 its
// trait value. Each species has a phylogenetic effect u (covariance Sigma_u (x) A) and a
// species effect v (covariance Sigma_v (x) I). The among-species standard deviation of each
// response, tau, is split by lambda, the phylogenetic share of the variance, which gets a
// beta prior directly: sd_u = tau * sqrt(lambda), sd_v = tau * sqrt(1 - lambda).
data {
  int<lower=1> S;                       // species
  matrix[S, S] L_A;                     // Cholesky factor of the phylogenetic correlation
  int<lower=0> N_p;                     // populations with a site mean
  vector[N_p] y_p;
  array[N_p] int<lower=1, upper=S> sp_p;
  vector<lower=0>[N_p] w_p;
  int<lower=1> N_t;                     // populations with a trait value
  vector[N_t] y_t;
  array[N_t] int<lower=1, upper=S> sp_t;
  vector<lower=0>[N_t] w_t;
  vector<lower=0>[2] lambda_prior;      // beta(a, b) for both lambdas
  real center;                          // standardisation of both responses
  real<lower=0> scale;
}
parameters {
  vector[2] mu;
  vector<lower=0>[2] tau;
  vector<lower=0, upper=1>[2] lambda;
  cholesky_factor_corr[2] L_u;
  cholesky_factor_corr[2] L_v;
  matrix[2, S] z_u;
  matrix[2, S] z_v;
  vector<lower=0>[2] sigma;
}
transformed parameters {
  vector[2] sd_u = tau .* sqrt(lambda);
  vector[2] sd_v = tau .* sqrt(1 - lambda);
  matrix[2, S] effect = rep_matrix(mu, S)
                        + diag_pre_multiply(sd_u, L_u) * z_u * L_A'
                        + diag_pre_multiply(sd_v, L_v) * z_v;
}
model {
  mu ~ normal(0, 1);
  tau ~ exponential(1);
  lambda ~ beta(lambda_prior[1], lambda_prior[2]);
  L_u ~ lkj_corr_cholesky(2);
  L_v ~ lkj_corr_cholesky(2);
  to_vector(z_u) ~ std_normal();
  to_vector(z_v) ~ std_normal();
  sigma ~ exponential(1);
  // weighted normal likelihood, as brms weights()
  if (N_p > 0) {
    vector[N_p] r = (y_p - effect[1, sp_p]') / sigma[1];
    target += -dot_product(w_p, 0.5 * square(r) + log(sigma[1]));
  }
  {
    vector[N_t] r = (y_t - effect[2, sp_t]') / sigma[2];
    target += -dot_product(w_t, 0.5 * square(r) + log(sigma[2]));
  }
}
generated quantities {
  real rho_u = L_u[2, 1];               // phylogenetic correlation of beta and the mean
  real rho_v = L_v[2, 1];
  real s2 = square(scale);
  real var_phylo_beta = square(sd_u[1]) * s2;
  real var_species_beta = square(sd_v[1]) * s2;
  real var_phylo_mean = square(sd_u[2]) * s2;
  real var_species_mean = square(sd_v[2]) * s2;
  real var_phylo_alpha = (square(sd_u[2]) + square(sd_u[1]) - 2 * rho_u * sd_u[1] * sd_u[2]) * s2;
  real var_species_alpha = (square(sd_v[2]) + square(sd_v[1]) - 2 * rho_v * sd_v[1] * sd_v[2]) * s2;
  real lambda_beta = lambda[1];
  real lambda_mean = lambda[2];
  real lambda_alpha = var_phylo_alpha / (var_phylo_alpha + var_species_alpha);
  real cor_phylo_beta_mean = rho_u;
  real cor_phylo_beta_alpha = (rho_u * sd_u[1] * sd_u[2] - square(sd_u[1])) * s2
                              / sqrt(var_phylo_beta * var_phylo_alpha);
  vector[S] beta = center + scale * effect[1]';
  vector[S] mean_trait = center + scale * effect[2]';
  vector[S] alpha = mean_trait - beta;
  vector[N_p] y_rep_p;
  vector[N_t] y_rep_t;
  vector[N_p] log_lik_p;
  vector[N_t] log_lik_t;
  for (i in 1:N_p) {
    real m = effect[1, sp_p[i]];
    y_rep_p[i] = center + scale * normal_rng(m, sigma[1]);
    log_lik_p[i] = w_p[i] * normal_lpdf(y_p[i] | m, sigma[1]);
  }
  for (i in 1:N_t) {
    real m = effect[2, sp_t[i]];
    y_rep_t[i] = center + scale * normal_rng(m, sigma[2]);
    log_lik_t[i] = w_t[i] * normal_lpdf(y_t[i] | m, sigma[2]);
  }
}
