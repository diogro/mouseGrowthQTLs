data {
  int<lower=1> K; // traits
  int<lower=1> J; // loci
  int<lower=0> N; // individuals
  int n_family;
  int family[N];
  vector[J] ad[N];
  vector[J] dm[N];
  vector[K] y[N];
}

transformed data {
  vector[K] zeros = rep_vector(0.0, K);
  real m0 = 10;           // Expected number of large slopes
  real slab_scale = 1;    // Scale for large slopes
  real slab_scale2 = square(slab_scale);
  real slab_df = 25;      // Effective degrees of freedom for large slopes
  real half_slab_df = 0.5 * slab_df;
}

parameters {
  // pre-shrinked slopes
  matrix[K, J] beta_tilde_ad;
  matrix[K, J] beta_tilde_dm;
  
  // global shrinkage scale, one per trait
  vector<lower=0>[K] tau_tilde_ad;
  vector<lower=0>[K] tau_tilde_dm;

  // local shrinkage scale, one per slope
  matrix<lower=0>[K, J] lambda_ad;
  matrix<lower=0>[K, J] lambda_dm;

  // good question, something like a threashold (maybe one per trait?!?)
  real<lower=0> c2_tilde_ad;
  real<lower=0> c2_tilde_dm;

  // intercept 
  vector[K] alpha;
  
  // Family means
  vector[K] beta_family[n_family];
  
  // G matrix
  cholesky_factor_corr[K] L_Omega_G;
  vector<lower=0>[K] L_sigma_G;
  
  // residual covariance
  cholesky_factor_corr[K] L_Omega_R; // cholesky correlation
  vector<lower=0>[K] L_sigma_R;      // standard deviations
}

transformed parameters {
  // individual means
  vector[K] mu[N];
  
  // global shrinkage scale, one per trait
  vector[K] tau0_ad;
  vector[K] tau_ad;
  vector[K] tau0_dm;
  vector[K] tau_dm;
  
  real c2_ad;
  real c2_dm;
  
  //local shrinkage scale
  matrix[K, J] lambda_tilde_ad;
  matrix[K, J] lambda_tilde_dm;
  
  // shrinked slopes
  matrix[K, J] beta_ad;
  matrix[K, J] beta_dm;
  {
    for(k in 1:K){
      tau0_ad[k] = (m0 / (J - m0)) * (L_sigma_R[k] / sqrt(1.0 * N));
      tau0_dm[k] = (m0 / (J - m0)) * (L_sigma_R[k] / sqrt(1.0 * N));
    }
    tau_ad = tau0_ad .* tau_tilde_ad; // tau_ad ~ cauchy(0, tau0_ad)
    tau_dm = tau0_dm .* tau_tilde_dm; // tau_ad ~ cauchy(0, tau0_ad)

    // c2_ad ~ inv_gamma(half_slab_df, half_slab_df * slab_scale2)
    // Implies that marginally beta_ad ~ student_t(slab_df, 0, slab_scale)
    c2_ad = slab_scale2 * c2_tilde_ad;
    c2_dm = slab_scale2 * c2_tilde_dm;

    for(j in 1:J){
      for(k in 1:K){
        lambda_tilde_ad[k, j] = tau_ad[k] * sqrt( c2_ad * square(lambda_ad[k, j]) ./ (c2_ad + square(tau_ad[k]) * square(lambda_ad[k, j])));
        lambda_tilde_dm[k, j] = tau_dm[k] * sqrt( c2_dm * square(lambda_dm[k, j]) ./ (c2_dm + square(tau_dm[k]) * square(lambda_dm[k, j])));
      }
    }
    // beta_ad ~ normal(0, tau_ad * lambda_tilde_ad)
    beta_ad = lambda_tilde_ad .* beta_tilde_ad;
    beta_dm = lambda_tilde_dm .* beta_tilde_dm;
  }
  
  for (n in 1:N)
    mu[n] = alpha + beta_ad * ad[n] + beta_dm * dm[n] + beta_family[family[n]];
}

model {
  matrix[K,K] L_Sigma_R;
  matrix[K,K] L_Sigma_G;
  
  L_Sigma_G = diag_pre_multiply(L_sigma_G, L_Omega_G);
  for (j in 1:n_family)
        beta_family[j] ~ multi_normal_cholesky(zeros, L_Sigma_G);

  L_Sigma_R = diag_pre_multiply(L_sigma_R, L_Omega_R);

  y ~ multi_normal_cholesky(mu, L_Sigma_R);
  
  c2_tilde_ad ~ inv_gamma(half_slab_df, half_slab_df);
  c2_tilde_dm ~ inv_gamma(half_slab_df, half_slab_df);
  
  to_vector(beta_tilde_ad) ~ normal(0, 1);
  to_vector(beta_tilde_dm) ~ normal(0, 1);
  to_vector(lambda_ad) ~ cauchy(0, 1);
  to_vector(lambda_dm) ~ cauchy(0, 1);
  tau_tilde_ad ~ cauchy(0, 1);
  tau_tilde_dm ~ cauchy(0, 1);
  
  alpha ~ normal(0, 2);
  L_Omega_R ~ lkj_corr_cholesky(4);
  L_sigma_R ~ normal(0, 1);
  
  L_Omega_G ~ lkj_corr_cholesky(4);
  L_sigma_G ~ normal(0, 1);
}

generated quantities {
    matrix[K, K] G;
    matrix[K, K] R;
    corr_matrix[K] corrG;
    corr_matrix[K] corrR;

    G = multiply_lower_tri_self_transpose(diag_pre_multiply(L_sigma_G, L_Omega_G));
    R = multiply_lower_tri_self_transpose(diag_pre_multiply(L_sigma_R, L_Omega_R));

    corrG = multiply_lower_tri_self_transpose(L_Omega_G);
    corrR = multiply_lower_tri_self_transpose(L_Omega_R);
}
