functions {
    // square root of a matrix (elementwise)
    matrix sqrt_mat(matrix x) {
        matrix[dims(x)[1], dims(x)[2]] res;
        for (n in 1:dims(x)[2]){
            for (m in 1:dims(x)[1]){
                res[m, n] = sqrt(x[m, n]);
            }
        }
        return res;
    }
    // square root of a vector (elementwise)
    vector sqrt_vec(vector x) {
        vector[dims(x)[1]] res;
        for (m in 1:dims(x)[1]){
            res[m] = sqrt(x[m]);
        }
        return res;
    }
}

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
  real slab_scale = 3;    // Scale for large slopes
  real slab_df = 25;      // Effective degrees of freedom for large slopes
  real nu_local = 1;     // degrees of freedom of lambda
  real scale_global = 1;
  real nu_global = 1;
}

parameters {
    // Horseshoe prior
    vector<lower=0>[K] r1_global_ad;
    vector<lower=0>[K] r1_global_dm;

    vector<lower=0>[K] r2_global_ad;
    vector<lower=0>[K] r2_global_dm;

    matrix<lower=0>[K, J] r1_local_ad;
    matrix<lower=0>[K, J] r2_local_ad;

    matrix<lower=0>[K, J] r1_local_dm;
    matrix<lower=0>[K, J] r2_local_dm;

    // what ever this is...
    real <lower=0> caux_ad;
    real <lower=0> caux_dm;

    // Effects matrix
    matrix[K,J] beta_ad;
    matrix[K,J] beta_dm;

    // Intercept
    vector[K] w0;

    // Family means
    vector[K] beta_family[n_family];

    // G matrix
    cholesky_factor_corr[K] L_Omega_G;
    vector<lower=0>[K] L_sigma_G;

    // R matrix
    cholesky_factor_corr[K] L_Omega_R;
    vector<lower=0>[K] r1_sigma_R;
    vector<lower=0>[K] r2_sigma_R;
}

transformed parameters{
    // global and local variance parameters, and the input weights
    vector<lower=0>[K] L_sigma_R;
    vector<lower=0>[K] tau_ad;
    vector<lower=0>[K] tau_dm;
    matrix<lower=0>[K, J] lambda_ad;
    matrix<lower=0>[K, J] lambda_dm;
    real < lower =0 > c_ad; // slab scale
    real < lower =0 > c_dm; // slab scale
    matrix<lower=0>[K, J] sd_theta_ad;
    matrix<lower=0>[K, J] sd_theta_dm;
    matrix[K, J] w_ad;
    matrix[K, J] w_dm;
    matrix[K,K] L_Sigma_R;
    vector[K] mu[N];

    L_sigma_R = 2.5 * r1_sigma_R .* sqrt_vec(r2_sigma_R);
    L_Sigma_R = diag_pre_multiply(L_sigma_R, L_Omega_R);

    tau_ad = r1_global_ad .* sqrt_vec(r2_global_ad) .* L_sigma_R * scale_global;
    tau_dm = r1_global_dm .* sqrt_vec(r2_global_dm) .* L_sigma_R * scale_global;

    lambda_ad = r1_local_ad .* sqrt_mat(r2_local_ad);
    lambda_dm = r1_local_dm .* sqrt_mat(r2_local_dm);

    c_ad = slab_scale * sqrt(caux_ad);
    c_dm = slab_scale * sqrt(caux_dm);
    
    for(j in 1:J){
        for(k in 1:K){
            sd_theta_ad[k, j] = tau_ad[k] * sqrt(c_ad^2 * square (lambda_ad[k, j]) ./ (c_ad^2 + tau_ad[k]^2 * square(lambda_ad[k, j])));
            sd_theta_dm[k, j] = tau_dm[k] * sqrt(c_dm^2 * square (lambda_dm[k, j]) ./ (c_dm^2 + tau_dm[k]^2 * square(lambda_dm[k, j])));
        }
    }
    w_ad = beta_ad .* sd_theta_ad;
    w_dm = beta_dm .* sd_theta_dm;

    for (n in 1:N)
        mu[n] = w0 + w_ad * ad[n] + w_dm * dm[n] + beta_family[family[n]];
}

model {
    matrix[K,K] L_Sigma_G;

    L_Sigma_G = diag_pre_multiply(L_sigma_G, L_Omega_G);

    for (j in 1:n_family)
        beta_family[j] ~ multi_normal_cholesky(zeros, L_Sigma_G);

    y ~ multi_normal_cholesky(mu, L_Sigma_R);

    //// half t-priors for lambdas (nu = 1 corresponds to horseshoe)
    to_vector(beta_ad) ~ normal(0, 1);
    to_vector(beta_dm) ~ normal(0, 1);

    to_vector(r1_local_ad) ~ normal(0.0, 1.0);
    to_vector(r2_local_ad) ~ inv_gamma(0.5*nu_local, 0.5*nu_local);

    to_vector(r1_local_dm) ~ normal(0.0, 1.0);
    to_vector(r2_local_dm) ~ inv_gamma(0.5*nu_local, 0.5*nu_local);
    // half cauchy for tau
    r1_global_ad ~ normal(0.0, 1.0);
    r2_global_ad ~ inv_gamma(0.5*nu_global, 0.5*nu_global);

    r1_global_dm ~ normal(0.0, 1.0);
    r2_global_dm ~ inv_gamma(0.5*nu_global, 0.5*nu_global);

    caux_ad ~ inv_gamma(0.5*slab_df , 0.5*slab_df);
    caux_dm ~ inv_gamma(0.5*slab_df , 0.5*slab_df);
  
    // weakly informative prior for the intercept
    w0 ~ normal(0,5);

    L_Omega_G ~ lkj_corr_cholesky(2);
    L_sigma_G ~ cauchy(0, 2.5);

    L_Omega_R ~ lkj_corr_cholesky(2);
    r1_sigma_R ~ normal(0.0, 1.0);
    r2_sigma_R ~ inv_gamma(0.5, 0.5);
}
generated quantities {
    matrix[K, K] G;
    matrix[K, K] R;
    corr_matrix[K] corrG;
    corr_matrix[K] corrR;
    matrix[K, J] shrink_ad;
    matrix[K, J] shrink_dm;

    G = multiply_lower_tri_self_transpose(diag_pre_multiply(L_sigma_G, L_Omega_G));
    R = multiply_lower_tri_self_transpose(diag_pre_multiply(L_sigma_R, L_Omega_R));

    corrG = multiply_lower_tri_self_transpose(L_Omega_G);
    corrR = multiply_lower_tri_self_transpose(L_Omega_R);

    for(j in 1:J){
        for(k in 1:K){
            shrink_ad[k, j] = 1 - 1/(1 + (lambda_ad[k, j]^2));
            shrink_dm[k, j] = 1 - 1/(1 + (lambda_dm[k, j]^2));
        }
    }
}

