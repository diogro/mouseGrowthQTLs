data {
    int<lower=1> K; // traits
    int<lower=0> N; // individuals
    int n_family;
    int family[N];
    vector[K] y[N];
}

transformed data{
    vector[K] zeros;
    zeros = rep_vector(0.0, K);
}

parameters {
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

model {
    vector[K] mu[N];
    matrix[K,K] L_Sigma_G;
    matrix[K,K] L_Sigma_R;

    L_Sigma_G = diag_pre_multiply(L_sigma_G, L_Omega_G);

    for (j in 1:n_family)
        beta_family[j] ~ multi_normal_cholesky(zeros, L_Sigma_G);

    L_Sigma_R = diag_pre_multiply(L_sigma_R, L_Omega_R);

    for (n in 1:N)
        mu[n] = w0 + beta_family[family[n]];

    y ~ multi_normal_cholesky(mu, L_Sigma_R);

    // weakly informative prior for the intercept
    w0 ~ normal(0,1);

    L_Omega_G ~ lkj_corr_cholesky(2);
    L_sigma_G ~ cauchy(0, 2.5);

    L_Omega_R ~ lkj_corr_cholesky(2);
    L_sigma_R ~ cauchy(0, 2.5);
}
generated quantities {
    matrix[K, K] G;
    matrix[K, K] R;
    matrix[K, K] P;
    corr_matrix[K] corrG;
    corr_matrix[K] corrR;

    G = multiply_lower_tri_self_transpose(diag_pre_multiply(L_sigma_G, L_Omega_G));
    R = multiply_lower_tri_self_transpose(diag_pre_multiply(L_sigma_R, L_Omega_R));
    P = G + R;

    corrG = multiply_lower_tri_self_transpose(L_Omega_G);
    corrR = multiply_lower_tri_self_transpose(L_Omega_R);
}
