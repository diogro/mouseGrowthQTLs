data{
    int N;
    int<lower=2> K;
    vector[K] mu;
    cov_matrix[K] Sigma;
    vector[K] growth[N];
    int FAMILY[N];
    int N_FAMILY;
    int addi[N];
    int domi[N];
}

transformed data{
    vector[K] zeros_FAMILY;
    matrix[K,K] L;

    L <- cholesky_decompose(Sigma);

    for ( i in 1:K ){
        zeros_FAMILY[i] <- 0;
    }
}

parameters{
    vector[K] alpha;
    vector[K] alpha_addi;
    vector[K] alpha_domi;
    vector[K] vary_FAMILY[N_FAMILY];
    cov_matrix[K] Sigma_FAMILY;
    cov_matrix[K] Sigma_Total;
}

transformed parameters {
    vector[K] Intercept;
    vector[K] beta_addi;
    vector[K] beta_domi;
    Intercept <- mu + L * alpha;
    beta_addi  <- mu + L * alpha_addi;
    beta_domi  <- mu + L * alpha_domi;
}

model{
    vector[K] vary[N];
    vector[K] glm[N];
    matrix[K, K] L_sigma_total;
    matrix[K, K] L_sigma_FAMILY;

    // Priors
    //Intercept ~ multi_normal( mu , Sigma );
    alpha ~ normal(0,1);
    alpha_addi ~ normal(0,1);
    alpha_domi ~ normal(0,1);

    L_sigma_FAMILY <- cholesky_decompose(Sigma_FAMILY);

    // Varying effects
    for ( j in 1:N_FAMILY ) vary_FAMILY[j] ~ multi_normal_cholesky( zeros_FAMILY , L_sigma_FAMILY );

    // Fixed effects

    L_sigma_total <- cholesky_decompose(Sigma_Total);

    for ( i in 1:N ) {
        vary[i] <- vary_FAMILY[FAMILY[i]];
        glm[i] <- vary[i] + Intercept + beta_addi * addi[i] + beta_domi * domi[i];
        growth[i] ~ multi_normal_cholesky( glm[i] , L_sigma_total );
    }
}
