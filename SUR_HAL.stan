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
    // power 2 of a vector of a vector over 2(elementwise)
    vector pow2vec(vector x) {
        vector[dims(x)[1]] res;
        for (m in 1:dims(x)[1]){
            res[m] = (x[m]*x[m])/2;
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
    # HAL prior
    vector<lower=0>[K] tau_ad;
    vector<lower=0>[K] tau_dm;


    matrix<lower=0>[K, J] lambda_ad;
    matrix<lower=0>[K, J] lambda_dm;

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
    matrix[K, J] w_ad;
    matrix[K, J] w_dm;
    vector<lower=0>[K] tau2_ad;
    vector<lower=0>[K] tau2_dm;

    tau2_ad = pow2vec(tau_ad);
    tau2_dm = pow2vec(tau_dm);

    w_ad = beta_ad .* lambda_ad;
    w_dm = beta_dm .* lambda_dm;
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

    to_vector(beta_ad) ~ normal(0, 1);
    to_vector(beta_dm) ~ normal(0, 1);

    tau_ad ~ inv_gamma(1, 1);
    tau_dm ~ inv_gamma(1, 1);

    for(k in 1:K){
        lambda_ad[k] ~ gamma(1, tau2_ad[k]);
        lambda_dm[k] ~ gamma(1, tau2_dm[k]);
    }

    // weakly informative prior for the intercept
    w0 ~ normal(0,5);

    L_Omega_G ~ lkj_corr_cholesky(2);
    L_sigma_G ~ cauchy(0, 2.5);

    L_Omega_R ~ lkj_corr_cholesky(2);
    L_sigma_R ~ cauchy(0, 2.5);
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
