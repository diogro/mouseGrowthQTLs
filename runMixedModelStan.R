setwd("/home/diogro/projects/mouse-qtls")
source('stanFunctions.R')
source('read_mouse_data.R')
source('utils.R')

library(MCMCglmm)

Rdatas_folder = "~/gdrive/LGSM_project_Rdatas/"
#Rdatas_folder = "./data/Rdatas/"

options(mc.cores = 1)

weight_data = inner_join(weight_phen_std, simulated_genomes[[8]], by = "ID")

weight_data = inner_join(weight_phen_std, markers, by = "ID")

weight_data[weight_traits] = scale(weight_data[weight_traits])

getStanInputMM = function(current_data, trait_vector)
{
    K        = length(trait_vector)
    N        = dim(current_data)[1]
    n_family = length(unique(current_data$FAMILY))
    family   = as.integer(as.factor(current_data$FAMILY))
    y        = as.matrix(current_data[trait_vector])
    param_list = list(K        = K,
                      N        = N,
                      n_family = n_family,
                      family   = family,
                      y        = y)
    return(param_list)
}

stan_MM = stan(file = "./mixedModelGmatrix.stan",
               data = getStanInputMM(weight_data, weight_traits),
               chain=3, iter = 400)

P = colMeans(rstan::extract(stan_MM, pars = c("P"))[[1]])
G = colMeans(rstan::extract(stan_MM, pars = c("G"))[[1]])
R = colMeans(rstan::extract(stan_MM, pars = c("R"))[[1]])
P  - cov(weight_data[weight_traits])

load(paste0(Rdatas_folder, "weight_MCMC_null_model.Rdata"))
summary(weight_MCMC_null_model)

G_mcmc = apply(array(weight_MCMC_null_model$VCV[,1:(num_weight_traits*num_weight_traits)], dim = c(1000, num_weight_traits, num_weight_traits)), 2:3, median)
R_mcmc = apply(array(weight_MCMC_null_model$VCV[,-c(1:(num_weight_traits*num_weight_traits))], dim = c(1000, num_weight_traits, num_weight_traits)), 2:3, median)

P_mcmc = G_mcmc + R_mcmc

data.frame(stan = cov2cor(P)[lower.tri(P)],
           MCMC = cov2cor(P_mcmc)[lower.tri(P)],
           P = cor(weight_data[weight_traits])[lower.tri(P)]) %>%
gather(variable, value, stan:MCMC) %>%
ggplot(aes(P, value, group = variable, color = variable)) + geom_point() + geom_abline()
