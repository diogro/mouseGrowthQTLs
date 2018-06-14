if(!require(rstan)){install.packages("rstan"); library(rstan)}
if(!require(plyr)){install.packages("plyr"); library(plyr)}
if(!require(tidyverse)){install.packages("tidyverse"); library(tidyverse)}
if(!require(mvtnorm)){install.packages("mvtnorm"); library(mvtnorm)}
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

N = 1500
J = 50
K = 2
m0 = 10
ad = cbind(matrix(rnorm(K*m0, 2, 0.1), K, m0), matrix(0, K, J-m0))
dm = cbind(matrix(rnorm(K*m0, 1, 0.1), K, m0), matrix(0, K, J-m0))

c_dark <- c("#8F2727")
c_dark_highlight <- c("#7C0000")
par(mfrow = c(1, 2), mar = c(4, 4, 0.5, 0.5))
hist(ad, main="", col=c_dark, border=c_dark_highlight,
     xlab="True Slopes", yaxt='n', ylim=c(0, 40), ylab="",
     breaks=2*(-100:100)/100)
hist(dm, main="", col=c_dark, border=c_dark_highlight,
     xlab="True Slopes", yaxt='n', ylim=c(0, 40), ylab="",
     breaks=2*(-100:100)/100)

X_ad = matrix(sample(c(1, 0, -1), N*J, replace = TRUE), N, J)
X_dm = abs(abs(X_ad) - 1) 
e = matrix(rnorm(N*K, 0, 0.5), N, K)
alpha = c(1, 1)

families = sample(1:200, N, replace = TRUE)
beta_family = rmvnorm(200, sigma = matrix(c(0.8, 0.5, 0.5, 0.8), 2, 2))

y = alpha + X_ad %*% t(ad) + X_dm %*% t(dm) + beta_family[families,] + e
coef(lm(y ~ X_ad + X_dm))

stan_data = list(N = N,
                 K = K,
                 J = J,
                 n_family = 200,
                 family = families,
                 ad = X_ad,
                 dm = X_dm,
                 y =  y)
fhs_fit = stan("SUR_fHS.stan", data = stan_data, iter = 200, chains = 2, 
               control=list(adapt_delta=0.99, max_treedepth=15))
betas_ad =  rstan::extract(fhs_fit, pars = "beta_ad")[[1]]
beta_ad = aaply(betas_ad, c(2, 3), mean)

estimates = rbind(data.frame(coef = 1:J, beta = beta_ad[1,], true = ad[1,], dim = 1),
                  data.frame(coef = 1:J, beta = beta_ad[2,], true = ad[2,], dim = 2)) %>%
  gather(trait, value, beta:true)

ggplot(estimates, aes(coef, value, group = trait, color = trait)) + geom_point() + facet_wrap(~dim)

Gs =  rstan::extract(fhs_fit, pars = "G")[[1]]
Rs =  rstan::extract(fhs_fit, pars = "R")[[1]]
aaply(Gs, c(2, 3), mean)
aaply(Rs, c(2, 3), mean)
cov(beta_family)
cov(e)

growth_data = inner_join(growth_phen_std, growth_markers, by = "ID")
growth_data[growth_traits] = scale(growth_data[growth_traits])

fhs_growth_fit = stan("SUR_fHS.stan", data = getStanInput(1, growth_data, growth_traits), 
                      iter = 60, chains = 1, 
                      control=list(adapt_delta=0.99, max_treedepth=15))
