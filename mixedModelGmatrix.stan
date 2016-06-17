data {
    int<lower=1> K;
    int<lower=1> J;
    int<lower=0> N;
    int n_family;
    int family[N];
    vector[J] ad[N];
    vector[J] dm[N];
    vector[K] y[N];
}

transformed data{
    vector[K] zeros;
    for ( i in 1:K )
        zeros[i] <- 0;
}

parameters {
    matrix[K,J] beta_ad;
    matrix[K,J] beta_dm;

    vector[K] beta_family[n_family];

    cholesky_factor_corr[K] L_Omega_G;
    vector<lower=0>[K] L_sigma_G;

    cholesky_factor_corr[K] L_Omega_R;
    vector<lower=0>[K] L_sigma_R;
}

model {
    vector[K] mu[N];
    matrix[K,K] L_Sigma_G;
    matrix[K,K] L_Sigma_R;

    L_Sigma_G <- diag_pre_multiply(L_sigma_G, L_Omega_G);
    L_Omega_G ~ lkj_corr_cholesky(4);
    L_sigma_G ~ cauchy(0, 2.5);

    for (j in 1:n_family)
        beta_family[j] ~ multi_normal_cholesky(zeros, L_Sigma_G);

    for (n in 1:N)
        mu[n] <- beta_ad * ad[n] + beta_dm * dm[n] + beta_family[family[n]];

    L_Sigma_R <- diag_pre_multiply(L_sigma_R, L_Omega_R);

    to_vector(beta_ad) ~ normal(0, 5);
    to_vector(beta_dm) ~ normal(0, 5);

    L_Omega_R ~ lkj_corr_cholesky(4);
    L_sigma_R ~ cauchy(0, 2.5);

    y ~ multi_normal_cholesky(mu, L_Sigma_R);
}
