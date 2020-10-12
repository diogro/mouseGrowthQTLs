setwd("/home/diogro/projects/mouse-qtls")
source('stanFunctions.R')
source('read_mouse_data.R')
source('utils.R')
old.par = par()

if(!require(MCMCglmm)){install.packages("MCMCglmm"); library(MCMCglmm)}

#Rdatas_folder = "~/gdrive/LGSM_project_Rdatas/"
Rdatas_folder = "./Rdatas/"

ncores = 8
registerDoMC(ncores)
options(mc.cores = 8)
setMKLthreads(ncores)

growth_data = growth_phen_std
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

growth_phen_growth_phen_sd_W = c(1, growth_phen_sd)
stan_MM_with_week9 = stan(file = "./mixedModelGmatrix.stan",
               data = getStanInputMM(growth_data, growth_traits_fitness),
               chain=4, iter = 2000, warmup = 500, control = list(max_treedepth = 11))
P_stan_w = colMeans(rstan::extract(stan_MM_with_week9, pars = c("P"))[[1]]) * outer(growth_phen_growth_phen_sd_W, growth_phen_growth_phen_sd_W)
G_stan_w = colMeans(rstan::extract(stan_MM_with_week9, pars = c("G"))[[1]]) * outer(growth_phen_growth_phen_sd_W, growth_phen_growth_phen_sd_W)
Gs_stan_w = aaply(rstan::extract(stan_MM_with_week9, pars = c("G"))[[1]], 1, '*', outer(growth_phen_growth_phen_sd_W, growth_phen_growth_phen_sd_W))
R_stan_w = colMeans(rstan::extract(stan_MM_with_week9, pars = c("R"))[[1]]) * outer(growth_phen_growth_phen_sd_W, growth_phen_growth_phen_sd_W)

stan_MM = stan(file = "./mixedModelGmatrix.stan",
               data = getStanInputMM(growth_data, growth_traits),
               chain=3, iter = 2000, warmup = 500, control = list(max_treedepth = 11))
P_stan = colMeans(rstan::extract(stan_MM, pars = c("P"))[[1]]) * outer(growth_phen_sd, growth_phen_sd)
G_stan = colMeans(rstan::extract(stan_MM, pars = c("G"))[[1]]) * outer(growth_phen_sd, growth_phen_sd)
Gs_stan = aaply(rstan::extract(stan_MM, pars = c("G"))[[1]], 1, '*', outer(growth_phen_sd, growth_phen_sd))
R_stan = colMeans(rstan::extract(stan_MM, pars = c("R"))[[1]]) * outer(growth_phen_sd, growth_phen_sd)
#load(file = "./Rdatas/growth_CovMatrices.Rdata")
P = cov(growth_phen_std[growth_traits]) * outer(growth_phen_sd, growth_phen_sd)

growth_phen_xf = growth_phen %>% 
  filter(xfostpair == "y") %>%
  mutate(NbyD = paste(Dam, NURSE, sep = "_")) %>%
  select(ID:NURSE, NbyD, everything())
formula = paste0("cbind(", paste(growth_traits, collapse = ", "), ")", "~ trait - 1" )
prior = list(R = list(V = diag(num_growth_traits), n = 0.002),
             G = list(G1 = list(V = diag(num_growth_traits) * 0.02, n = 0.001),
                      G2 = list(V = diag(num_growth_traits) * 0.02, n = 0.001)))
growth_MCMC_null_model = MCMCglmm(as.formula(formula), 
                                  data = as.data.frame(growth_phen_xf), 
                                  random = ~us(trait):Dam + us(trait):NURSE,
                                  rcov = ~us(trait):ID,
                                  prior = prior,
                                  nitt = 103000, burnin = 3000, thin = 100,
                                  family = rep("gaussian", length(growth_traits)))

save(growth_MCMC_null_model, file = paste0(Rdatas_folder, "growth_MCMC_DamNurse_model.Rdata"))
load(paste0(Rdatas_folder, "growth_MCMC_DamNurse_model.Rdata"))
load(paste0(Rdatas_folder, "growth_MCMC_null_model.Rdata"))
summary(growth_MCMC_null_model)

id = colnames(growth_MCMC_null_model$VCV)
Gs_dam = array(growth_MCMC_null_model$VCV[,grep("Dam", id)], 
                    dim = c(1000, num_growth_traits, num_growth_traits))
G_dam = apply(array(growth_MCMC_null_model$VCV[,grep("Dam", id)], 
                     dim = c(1000, num_growth_traits, num_growth_traits)), 2:3, median)
G_nurse = apply(array(growth_MCMC_null_model$VCV[,grep("NURSE", id)], 
                     dim = c(1000, num_growth_traits, num_growth_traits)), 2:3, median)
#G_NbyD = apply(array(growth_MCMC_null_model$VCV[,grep("NbyD", id)], 
#                      dim = c(1000, num_growth_traits, num_growth_traits)), 2:3, median)
R_mcmc = apply(array(growth_MCMC_null_model$VCV[,grep("ID", id)], 
                     dim = c(1000, num_growth_traits, num_growth_traits)), 2:3, median)

png("./data/growth_fullSib_Dam_Nurse_Gcorrelation.png", width = 2100, height = 700)
par(mfrow = c(1, 3), cex=2, oma = c(0, 0, 0, 0))
corrplot.mixed(cov2cor(G_stan), upper = "ellipse", mar = c(0, 0, 1, 0), main = "FullSib")
corrplot.mixed(cov2cor(G_dam),   upper = "ellipse", mar = c(0, 0, 1, 0), main = "Dam")
corrplot.mixed(cov2cor(G_nurse), upper = "ellipse", mar = c(0, 0, 1, 0), main = "Nurse")
dev.off()

P_mcmc = G_dam + G_nurse + R_mcmc
G_mcmc = G_dam + G_nurse

corrplot.mixed(cov2cor(G_stan), upper = "ellipse")

data.frame(stan = (P_stan)[lower.tri(P, diag = TRUE)],
           MCMC = (P_mcmc)[lower.tri(P, diag = TRUE)],
           P = (P)[lower.tri(P, diag = TRUE)]) %>%
gather(variable, value, stan:MCMC) %>%
ggplot(aes(P, value, group = variable, color = variable)) + geom_point() + geom_abline()

## Cross-fosterd vs non-cross-fostered

growth_data_cf = filter(growth_data, Dam != NURSE)
growth_data_ncf = filter(growth_data, Dam == NURSE)

stan_cf = stan(file = "./mixedModelGmatrix.stan",
               data = getStanInputMM(growth_data_cf, growth_traits),
               chain=3, iter = 2000, warmup = 500, control = list(max_treedepth = 11, adapt_delta = 0.99))
stan_ncf = stan(file = "./mixedModelGmatrix.stan",
                data = getStanInputMM(growth_data_ncf, growth_traits),
                chain=3, iter = 2000, control = list(max_treedepth = 11, adapt_delta = 0.99))

G_cf  = colMeans(rstan::extract(stan_cf, pars = c("G"))[[1]]) * outer(growth_phen_sd, growth_phen_sd)
G_ncf = colMeans(rstan::extract(stan_ncf, pars = c("G"))[[1]]) * outer(growth_phen_sd, growth_phen_sd)

png("./data/growth_fullSib_CF_nCF_Gcorrelation.png", width = 2100, height = 700)
par(mfrow = c(1, 3), cex=2, oma = c(0, 0, 0, 0))
corrplot.mixed(cov2cor(G_cf),   upper = "ellipse", mar = c(0, 0, 1, 0), main = "FullSib CF")
corrplot.mixed(cov2cor(G_stan), upper = "ellipse", mar = c(0, 0, 1, 0), main = "FullSib Family")
corrplot.mixed(cov2cor(G_ncf),  upper = "ellipse", mar = c(0, 0, 1, 0), main = "FullSib Non-CF")
par(old.par)
dev.off()

save(P_stan_w, G_stan_w, Gs_stan_w, R_stan_w, 
     P_stan, G_stan, Gs_stan, R_stan, G_cf, G_ncf, file = "./Rdatas/growth_CovMatrices.Rdata")
load(file = "./Rdatas/growth_CovMatrices.Rdata")

### Dam nurse separation
# 
# growth_data = inner_join(growth_phen, growth_markers, by = "ID")
# getStanInputDamNurseMM = function(current_data, trait_vector)
# {
#   K        = length(trait_vector)
#   N        = dim(current_data)[1]
#   n_dam = length(unique(current_data$Dam))
#   dam   = as.integer(as.factor(current_data$Dam))
#   n_nurse = length(unique(current_data$NURSE))
#   nurse   = as.integer(as.factor(current_data$NURSE))
#   y        = as.matrix(current_data[trait_vector])
#   param_list = list(K       = K,
#                     N       = N,
#                     n_dam   = n_dam,
#                     dam     = dam,
#                     n_nurse = n_nurse,
#                     nurse   = nurse,
#                     y       = y)
#   return(param_list)
# }
# 
# stan_DamNurse = stan(file = "./mixedModelDamNurseGmatrix.stan",
#                data = getStanInputDamNurseMM(growth_data, growth_traits),
#                chain=8, iter = 2000, control = list(max_treedepth = 11, adapt_delta = 0.99))
# save(stan_DamNurse, file = "./Rdatas/growth_DamNurseGMatrixFit.Rdata")
# 
# 
# P_DM = colMeans(rstan::extract(stan_DamNurse, pars = c("P"))[[1]])
# G_dam = colMeans(rstan::extract(stan_DamNurse, pars = c("G_dam"))[[1]])
# G_nurse = colMeans(rstan::extract(stan_DamNurse, pars = c("G_nurse"))[[1]])
# Gs_dam = rstan::extract(stan_DamNurse, pars = c("G_dam"))[[1]]
# Gs_nurse = rstan::extract(stan_DamNurse, pars = c("G_nurse"))[[1]]
# R_DM = colMeans(rstan::extract(stan_DamNurse, pars = c("R"))[[1]])
# save(P_DM, G_dam, G_nurse, Gs_dam, Gs_nurse, R_DM, file = "./Rdatas/growth_DamNurseMatrices.Rdata")
# load(file = "./Rdatas/growth_DamNurseMatrices.Rdata")
# 
# 
# par(mfrow = c(1, 3), cex=2, oma = c(0, 0, 0, 0))
# corrplot.mixed(cov2cor(G_dam),   upper = "ellipse", mar = c(0, 0, 1, 0), main = "FullSib Dam")
# corrplot.mixed(cov2cor(G_stan), upper = "ellipse", mar = c(0, 0, 1, 0), main = "FullSib Family")
# corrplot.mixed(cov2cor(G_nurse), upper = "ellipse", mar = c(0, 0, 1, 0), main = "FullSib Nurse")



