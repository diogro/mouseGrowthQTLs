data {
    int<lower=1> K; // traits
    int<lower=0> N; // individuals
    int n_dam;
    int dam[N];
    int n_nurse;
    int nurse[N];
    vector[K] y[N];
}

transformed data{
    vector[K] zeros;
    zeros = rep_vector(0.0, K);
}

parameters {
    // Intercept
    vector[K] w0;

    // Dam means
    vector[K] beta_dam[n_dam];
    
    // Nurse means
    vector[K] beta_nurse[n_nurse];

    // Dam matrix
    cholesky_factor_corr[K] L_Omega_dam;
    vector<lower=0>[K] L_sigma_dam;

    // Nurse matrix
    cholesky_factor_corr[K] L_Omega_nurse;
    vector<lower=0>[K] L_sigma_nurse;

    // R matrix
    cholesky_factor_corr[K] L_Omega_R;
    vector<lower=0>[K] L_sigma_R;
}

model {
    vector[K] mu[N];
    matrix[K,K] L_Sigma_dam;
    matrix[K,K] L_Sigma_nurse;  
    matrix[K,K] L_Sigma_R;

    L_Sigma_dam = diag_pre_multiply(L_sigma_dam, L_Omega_dam);
    L_Sigma_nurse = diag_pre_multiply(L_sigma_nurse, L_Omega_nurse);

    for (j in 1:n_dam)
        beta_dam[j] ~ multi_normal_cholesky(zeros, L_Sigma_dam);
        
    for (j in 1:n_nurse)
        beta_nurse[j] ~ multi_normal_cholesky(zeros, L_Sigma_nurse);

    L_Sigma_R = diag_pre_multiply(L_sigma_R, L_Omega_R);

    for (n in 1:N)
        mu[n] = w0 + beta_dam[dam[n]] + beta_nurse[nurse[n]];

    y ~ multi_normal_cholesky(mu, L_Sigma_R);

    // weakly informative prior for the intercept
    w0 ~ normal(0,1);

    L_Omega_dam ~ lkj_corr_cholesky(2);
    L_sigma_dam ~ cauchy(0, 2.5);
    
    L_Omega_nurse ~ lkj_corr_cholesky(2);
    L_sigma_nurse ~ cauchy(0, 2.5);

    L_Omega_R ~ lkj_corr_cholesky(2);
    L_sigma_R ~ cauchy(0, 2.5);
}
generated quantities {
    matrix[K, K] G_dam;
    matrix[K, K] G_nurse;
    matrix[K, K] R;
    matrix[K, K] P;
    corr_matrix[K] corrG_dam;
    corr_matrix[K] corrG_nurse;
    corr_matrix[K] corrR;

    G_dam = multiply_lower_tri_self_transpose(diag_pre_multiply(L_sigma_dam, L_Omega_dam));
    G_nurse = multiply_lower_tri_self_transpose(diag_pre_multiply(L_sigma_nurse, L_Omega_nurse));
    R = multiply_lower_tri_self_transpose(diag_pre_multiply(L_sigma_R, L_Omega_R));
    P = G_dam + G_nurse + R;

    corrG_dam = multiply_lower_tri_self_transpose(L_Omega_dam);
    corrG_nurse = multiply_lower_tri_self_transpose(L_Omega_nurse);
    corrR = multiply_lower_tri_self_transpose(L_Omega_R);
}
