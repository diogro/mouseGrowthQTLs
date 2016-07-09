functions {
    // square root of a matrix (elementwise)
    matrix sqrt_vec(matrix x) {
        matrix[dims(x)[1], dims(x)[2]] res;
        for (n in 1:dims(x)[2]){
            for (m in 1:dims(x)[1]){
                res[m, n] = sqrt(x[m, n]);
            }
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

transformed data{
    vector[K] zeros;
    zeros = rep_vector(0.0, K);
}

parameters {
    # Horseshoe prior
    real<lower=0> r1_global;
    real<lower=0> r2_global;

    matrix<lower=0>[K, J] r1_local_ad;
    matrix<lower=0>[K, J] r2_local_ad;

    matrix<lower=0>[K, J] r1_local_dm;
    matrix<lower=0>[K, J] r2_local_dm;

    matrix<lower=0>[K, J] r1_localPlus_ad;
    matrix<lower=0>[K, J] r2_localPlus_ad;

    matrix<lower=0>[K, J] r1_localPlus_dm;
    matrix<lower=0>[K, J] r2_localPlus_dm;

    # Effects matrix
    matrix[K,J] beta_ad;
    matrix[K,J] beta_dm;

    # Intercept
    vector[K] w0;

    # Family means
    vector[K] beta_family[n_family];

    # G matrix
    cholesky_factor_corr[K] L_Omega_G;
    vector<lower=0>[K] L_sigma_G;

    # R matrix
    cholesky_factor_corr[K] L_Omega_R;
    vector<lower=0>[K] L_sigma_R;
}

transformed parameters{
    // global and local variance parameters, and the input weights
    real<lower=0> tau;
    matrix<lower=0>[K, J] lambda_ad;
    matrix<lower=0>[K, J] lambda_dm;
    matrix<lower=0>[K, J] lambdaPlus_ad;
    matrix<lower=0>[K, J] lambdaPlus_dm;
    matrix[K, J] w_ad;
    matrix[K, J] w_dm;

    tau = r1_global * sqrt(r2_global);

    lambda_ad = r1_local_ad .* sqrt_vec(r2_local_ad);
    lambda_dm = r1_local_dm .* sqrt_vec(r2_local_dm);

    lambdaPlus_ad = r1_localPlus_ad .* sqrt_vec(r2_localPlus_ad);
    lambdaPlus_dm = r1_localPlus_dm .* sqrt_vec(r2_localPlus_dm);

    w_ad = beta_ad .* lambda_ad .* lambdaPlus_ad * tau;
    w_dm = beta_dm .* lambda_dm .* lambdaPlus_dm * tau;
}

model {
    vector[K] mu[N];
    matrix[K,K] L_Sigma_G;
    matrix[K,K] L_Sigma_R;

    L_Sigma_G = diag_pre_multiply(L_sigma_G, L_Omega_G);

    for (j in 1:n_family)
        beta_family[j] ~ multi_normal_cholesky(zeros, L_Sigma_G);

    L_Sigma_R = diag_pre_multiply(L_sigma_R, L_Omega_R);

    for (n in 1:N)
        mu[n] = w0 + w_ad * ad[n] + w_dm * dm[n] + beta_family[family[n]];

    y ~ multi_normal_cholesky(mu, L_Sigma_R);

    #// half t-priors for lambdas (nu = 1 corresponds to horseshoe)
    to_vector(beta_ad) ~ normal(0, 1);
    to_vector(r1_local_ad) ~ normal(0.0, 1.0);
    to_vector(r2_local_ad) ~ inv_gamma(0.5*3, 0.5*3);
    to_vector(beta_dm) ~ normal(0, 1);
    to_vector(r1_local_dm) ~ normal(0.0, 1.0);
    to_vector(r2_local_dm) ~ inv_gamma(0.5*3, 0.5*3);
    to_vector(r1_localPlus_dm) ~ normal(0.0, 1.0);
    to_vector(r2_localPlus_dm) ~ inv_gamma(0.5*3, 0.5*3);
    // half cauchy for tau
    r1_global ~ normal(0.0, 1.0);
    r2_global ~ inv_gamma(0.5, 0.5);

    // weakly informative prior for the intercept 
    w0 ~ normal(0,5);

    L_Omega_G ~ lkj_corr_cholesky(4);
    L_sigma_G ~ cauchy(0, 2.5);

    L_Omega_R ~ lkj_corr_cholesky(4);
    L_sigma_R ~ cauchy(0, 2.5);
}
