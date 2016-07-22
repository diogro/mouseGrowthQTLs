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

transformed data{
    vector[K] zeros;
    zeros = rep_vector(0.0, K);
}

parameters {
    matrix[K, J] r1_local_ad;
    matrix<lower=0>[K, J] r2_local_ad;

    matrix[K, J] r1_local_dm;
    matrix<lower=0>[K, J] r2_local_dm;

    # Intercept
    vector[K] w0;

    # Family means
    vector[K] beta_family[n_family];

    # G matrix
    cholesky_factor_corr[K] L_Omega_G;
    vector<lower=0>[K] L_sigma_G;

    # R matrix
    cholesky_factor_corr[K] L_Omega_R;
    vector<lower=0>[K] r1_sigma_R;
    vector<lower=0>[K] r2_sigma_R;
}

transformed parameters{
    // global and local variance parameters, and the input weights
    vector<lower=0>[K] L_sigma_R;
    matrix[K, J] w_ad;
    matrix[K, J] w_dm;

    L_sigma_R = 2.5 * r1_sigma_R .* sqrt_vec(r2_sigma_R);

    w_ad = r1_local_ad .* sqrt_mat(r2_local_ad);
    w_dm = r1_local_dm .* sqrt_mat(r2_local_dm);
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

    to_vector(r1_local_ad) ~ normal(0.0, 1.0);
    to_vector(r2_local_ad) ~ inv_gamma(0.5*3, 0.5*3);
    to_vector(r1_local_dm) ~ normal(0.0, 1.0);
    to_vector(r2_local_dm) ~ inv_gamma(0.5*3, 0.5*3);

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

    G = multiply_lower_tri_self_transpose(diag_pre_multiply(L_sigma_G, L_Omega_G));
    R = multiply_lower_tri_self_transpose(diag_pre_multiply(L_sigma_R, L_Omega_R));

    corrG = multiply_lower_tri_self_transpose(L_Omega_G);
    corrR = multiply_lower_tri_self_transpose(L_Omega_R);
}

