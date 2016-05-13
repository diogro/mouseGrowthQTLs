data{
    int N;
    int N_locus;
    int<lower=2> K;
    row_vector[K] mu;
    cov_matrix[K] Sigma;
    vector[K] growth[N];
    int FAMILY[N];
    int N_FAMILY;
    matrix[N, N_locus] addi;
}

transformed data{
    vector[K] zeros_FAMILY;
    matrix[K,K] L;
    int max_QTL;

    max_QTL <- 1080;

    L <- cholesky_decompose(Sigma);

    for ( i in 1:K ){
        zeros_FAMILY[i] <- 0;
    }
}

parameters{
    row_vector[K] alpha;
    matrix[N_locus, K] beta_addi;
    vector[K] vary_FAMILY[N_FAMILY];
    cov_matrix[K] Sigma_FAMILY;
    cov_matrix[K] Sigma_Total;
    real<lower=0, upper=1> omega;
    real<lower=0, upper=max_QTL> alpha_lambda;
    real<lower=0, upper=max_QTL> beta_lambda;
}

transformed parameters {
    row_vector[K] Intercept;
    Intercept <- mu + alpha * L;
}

model{
    matrix[N, K] vary;
    matrix[N, K] locus;
    matrix[K, K] L_sigma_FAMILY;
    matrix[K, K] L_sigma_Total;
    vector[K] glm[N];
    int lambda[N_locus];

    // Indicator vector
    omega ~ beta(alpha_lambda,beta_lambda);
    for (n in 1:N_locus)
        lambda[n] ~ bernoulli(omega);

    // Priors
    //Intercept ~ multi_normal( mu , Sigma );
    alpha ~ normal(0,1);

    L_sigma_FAMILY <- cholesky_decompose(Sigma_FAMILY);

    // Varying effects
    for ( j in 1:N_FAMILY ) vary_FAMILY[j] ~ multi_normal_cholesky( zeros_FAMILY , L_sigma_FAMILY );

    // Fixed effects

    L_sigma_Total <- cholesky_decompose(Sigma_Total);

    locus <- addi * diag_matrix(lambda) * beta_addi
    for ( i in 1:N ) {
        vary[i] <- Intercept + vary_FAMILY[FAMILY[i]]';
        for(j in 1:K)
            glm[i,j] <- vary[i,j] + locus[i,j];
        growth[i] ~ multi_normal_cholesky( glm[i] , L_sigma_Total );
    }
}
