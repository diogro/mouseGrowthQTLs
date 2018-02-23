setwd("/home/diogro/projects/mouse-qtls")
source('stanFunctions.R')
source('read_mouse_data.R')
source('utils.R')

library(MCMCglmm)

#Rdatas_folder = "~/gdrive/LGSM_project_Rdatas/"
Rdatas_folder = "./Rdatas/"

options(mc.cores = 8)

growth_data = inner_join(growth_phen, growth_markers, by = "ID")
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
               data = getStanInputMM(growth_data, growth_traits),
               chain=6, iter = 600, control = list(max_treedepth = 11))

P_stan = colMeans(rstan::extract(stan_MM, pars = c("P"))[[1]])
G_stan = colMeans(rstan::extract(stan_MM, pars = c("G"))[[1]])
R_stan = colMeans(rstan::extract(stan_MM, pars = c("R"))[[1]])
save(P_stan, G_stan, R_stan, file = "./Rdatas/growth_CovMatrices.Rdata")
P = cov(growth_data[growth_traits])
cov2cor(G)

load(paste0(Rdatas_folder, "growth_MCMC_null_model.Rdata"))
summary(growth_MCMC_null_model)

G_mcmc = apply(array(growth_MCMC_null_model$VCV[,1:(num_growth_traits*num_growth_traits)], dim = c(1000, num_growth_traits, num_growth_traits)), 2:3, median)
R_mcmc = apply(array(growth_MCMC_null_model$VCV[,-c(1:(num_growth_traits*num_growth_traits))], dim = c(1000, num_growth_traits, num_growth_traits)), 2:3, median)

P_mcmc = G_mcmc + R_mcmc

library(corrplot)
corrplot.mixed(cov2cor(G_stan), upper = "ellipse")

data.frame(stan = (P_stan)[lower.tri(P)],
           MCMC = (P_mcmc)[lower.tri(P)],
           P = cov(growth_data[growth_traits])[lower.tri(P)]) %>%
gather(variable, value, stan:MCMC) %>%
ggplot(aes(P, value, group = variable, color = variable)) + geom_point() + geom_abline()
