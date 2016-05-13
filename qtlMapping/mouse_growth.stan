data{
    int N;
    real value[N];
    real variablegrow23[N];
    real variablegrow34[N];
    real variablegrow45[N];
    real variablegrow56[N];
    real variablegrow67[N];
    real variablegrow78[N];
    real A1[N];
    real D1[N];
    int FAMILY[N];
    real variablegrow23_X_A1[N];
    real variablegrow34_X_A1[N];
    real variablegrow45_X_A1[N];
    real variablegrow56_X_A1[N];
    real variablegrow67_X_A1[N];
    real variablegrow78_X_A1[N];
    real variablegrow23_X_D1[N];
    real variablegrow34_X_D1[N];
    real variablegrow45_X_D1[N];
    real variablegrow56_X_D1[N];
    real variablegrow67_X_D1[N];
    real variablegrow78_X_D1[N];
    int N_FAMILY;
}

transformed data{
    vector[7] zeros_FAMILY;
    for ( i in 1:7 ) zeros_FAMILY[i] <- 0;
}

parameters{
    real Intercept;
    real beta_A1;
    real beta_D1;
    real beta_variablegrow23_X_A1;
    real beta_variablegrow34_X_A1;
    real beta_variablegrow45_X_A1;
    real beta_variablegrow56_X_A1;
    real beta_variablegrow67_X_A1;
    real beta_variablegrow78_X_A1;
    real beta_variablegrow23_X_D1;
    real beta_variablegrow34_X_D1;
    real beta_variablegrow45_X_D1;
    real beta_variablegrow56_X_D1;
    real beta_variablegrow67_X_D1;
    real beta_variablegrow78_X_D1;
    real<lower=0> sigma;
    vector[7] vary_FAMILY[N_FAMILY];
    cov_matrix[7] Sigma_FAMILY;
}

model{
    real vary[N];
    real glm[N];
    // Priors
    Intercept ~ normal( 0 , 100 );
    beta_A1 ~ normal( 0 , 100 );
    beta_D1 ~ normal( 0 , 100 );
    beta_variablegrow23_X_A1 ~ normal( 0 , 100 );
    beta_variablegrow34_X_A1 ~ normal( 0 , 100 );
    beta_variablegrow45_X_A1 ~ normal( 0 , 100 );
    beta_variablegrow56_X_A1 ~ normal( 0 , 100 );
    beta_variablegrow67_X_A1 ~ normal( 0 , 100 );
    beta_variablegrow78_X_A1 ~ normal( 0 , 100 );
    beta_variablegrow23_X_D1 ~ normal( 0 , 100 );
    beta_variablegrow34_X_D1 ~ normal( 0 , 100 );
    beta_variablegrow45_X_D1 ~ normal( 0 , 100 );
    beta_variablegrow56_X_D1 ~ normal( 0 , 100 );
    beta_variablegrow67_X_D1 ~ normal( 0 , 100 );
    beta_variablegrow78_X_D1 ~ normal( 0 , 100 );
    sigma ~ uniform( 0 , 100 );
    // Varying effects
    for ( j in 1:N_FAMILY ) vary_FAMILY[j] ~ multi_normal( zeros_FAMILY , Sigma_FAMILY );
    // Fixed effects
    for ( i in 1:N ) {
        vary[i] <- vary_FAMILY[FAMILY[i],1]
                + vary_FAMILY[FAMILY[i],2] * variablegrow23[i]
                + vary_FAMILY[FAMILY[i],3] * variablegrow34[i]
                + vary_FAMILY[FAMILY[i],4] * variablegrow45[i]
                + vary_FAMILY[FAMILY[i],5] * variablegrow56[i]
                + vary_FAMILY[FAMILY[i],6] * variablegrow67[i]
                + vary_FAMILY[FAMILY[i],7] * variablegrow78[i];
        glm[i] <- vary[i] + Intercept
                + beta_A1 * A1[i]
                + beta_D1 * D1[i]
                + beta_variablegrow23_X_A1 * variablegrow23_X_A1[i]
                + beta_variablegrow34_X_A1 * variablegrow34_X_A1[i]
                + beta_variablegrow45_X_A1 * variablegrow45_X_A1[i]
                + beta_variablegrow56_X_A1 * variablegrow56_X_A1[i]
                + beta_variablegrow67_X_A1 * variablegrow67_X_A1[i]
                + beta_variablegrow78_X_A1 * variablegrow78_X_A1[i]
                + beta_variablegrow23_X_D1 * variablegrow23_X_D1[i]
                + beta_variablegrow34_X_D1 * variablegrow34_X_D1[i]
                + beta_variablegrow45_X_D1 * variablegrow45_X_D1[i]
                + beta_variablegrow56_X_D1 * variablegrow56_X_D1[i]
                + beta_variablegrow67_X_D1 * variablegrow67_X_D1[i]
                + beta_variablegrow78_X_D1 * variablegrow78_X_D1[i];
    }
    value ~ normal( glm , sigma );
}
