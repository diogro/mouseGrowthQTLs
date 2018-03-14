setwd("/home/diogro/projects/mouse-qtls")
source('stanFunctions.R')
source('read_mouse_data.R')
source('utils.R')
old.par = par()

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
Gs_stan = rstan::extract(stan_MM, pars = c("G"))[[1]]
R_stan = colMeans(rstan::extract(stan_MM, pars = c("R"))[[1]])
save(P_stan, G_stan, Gs_stan, R_stan, file = "./Rdatas/growth_CovMatrices.Rdata")
load(file = "./Rdatas/growth_CovMatrices.Rdata")
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

### Dam nurse separation

growth_data = inner_join(growth_phen, growth_markers, by = "ID")
getStanInputDamNurseMM = function(current_data, trait_vector)
{
  K        = length(trait_vector)
  N        = dim(current_data)[1]
  n_dam = length(unique(current_data$Dam))
  dam   = as.integer(as.factor(current_data$Dam))
  n_nurse = length(unique(current_data$NURSE))
  nurse   = as.integer(as.factor(current_data$NURSE))
  y        = as.matrix(current_data[trait_vector])
  param_list = list(K       = K,
                    N       = N,
                    n_dam   = n_dam,
                    dam     = dam,
                    n_nurse = n_nurse,
                    nurse   = nurse,
                    y       = y)
  return(param_list)
}

stan_DamNurse = stan(file = "./mixedModelDamNurseGmatrix.stan",
               data = getStanInputDamNurseMM(growth_data, growth_traits),
               chain=8, iter = 2000, control = list(max_treedepth = 11, adapt_delta = 0.99))
save(stan_DamNurse, file = "./Rdatas/growth_DamNurseGMatrixFit.Rdata")


P_DM = colMeans(rstan::extract(stan_DamNurse, pars = c("P"))[[1]])
G_dam = colMeans(rstan::extract(stan_DamNurse, pars = c("G_dam"))[[1]])
G_nurse = colMeans(rstan::extract(stan_DamNurse, pars = c("G_nurse"))[[1]])
Gs_dam = rstan::extract(stan_DamNurse, pars = c("G_dam"))[[1]]
Gs_nurse = rstan::extract(stan_DamNurse, pars = c("G_nurse"))[[1]]
R_DM = colMeans(rstan::extract(stan_DamNurse, pars = c("R"))[[1]])
save(P_DM, G_dam, G_nurse, Gs_dam, Gs_nurse, R_DM, file = "./Rdatas/growth_DamNurseMatrices.Rdata")
load(file = "./Rdatas/growth_DamNurseMatrices.Rdata")

png("./data/growth_fullSib_Dam_Nurse_Gcorrelation.png", width = 2100, height = 700)
par(mfrow = c(1, 3), cex=2, oma = c(0, 0, 0, 0))
corrplot.mixed(cov2cor(G_dam),   upper = "ellipse", mar = c(0, 0, 1, 0), main = "FullSib Dam")
corrplot.mixed(cov2cor(G_stan), upper = "ellipse", mar = c(0, 0, 1, 0), main = "FullSib Family")
corrplot.mixed(cov2cor(G_nurse), upper = "ellipse", mar = c(0, 0, 1, 0), main = "FullSib Nurse")
par(old.par)
dev.off()

