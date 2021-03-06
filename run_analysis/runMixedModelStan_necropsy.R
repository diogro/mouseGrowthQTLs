setwd("/home/diogro/projects/mouse-qtls")
source('stanFunctions.R')
source('read_mouse_data.R')
source('utils.R')
old.par = par()

if(!require(MCMCglmm)){install.packages("MCMCglmm"); library(MCMCglmm)}

#Rdatas_folder = "~/gdrive/LGSM_project_Rdatas/"
Rdatas_folder = "./Rdatas/"

ncores = 4
registerDoMC(ncores)
options(mc.cores = ncores)
setMKLthreads(ncores)

necropsy_data = necropsy_phen_std
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
               data = getStanInputMM(necropsy_data, necropsy_traits),
               chain=4, iter = 2000, warmup = 500, control = list(max_treedepth = 11))
P_stan = colMeans(rstan::extract(stan_MM, pars = c("P"))[[1]]) * outer(necropsy_phen_sd, necropsy_phen_sd)
G_stan = colMeans(rstan::extract(stan_MM, pars = c("G"))[[1]]) * outer(necropsy_phen_sd, necropsy_phen_sd)
Gs_stan = aaply(rstan::extract(stan_MM, pars = c("G"))[[1]], 1, '*', outer(necropsy_phen_sd, necropsy_phen_sd))
R_stan = colMeans(rstan::extract(stan_MM, pars = c("R"))[[1]]) * outer(necropsy_phen_sd, necropsy_phen_sd)
#load(file = "./Rdatas/necropsy_CovMatrices.Rdata")
P = cov(necropsy_phen_std[necropsy_traits]) * outer(necropsy_phen_sd, necropsy_phen_sd)

necropsy_phen_xf = necropsy_phen %>% 
  filter(xfostpair == "y") %>%
  mutate(NbyD = paste(Dam, NURSE, sep = "_")) %>%
  select(ID:NURSE, NbyD, everything())
formula = paste0("cbind(", paste(necropsy_traits, collapse = ", "), ")", "~ trait - 1" )
prior = list(R = list(V = diag(num_necropsy_traits), n = 0.002),
             G = list(G1 = list(V = diag(num_necropsy_traits) * 0.02, n = 0.001),
                      G2 = list(V = diag(num_necropsy_traits) * 0.02, n = 0.001)))
necropsy_MCMC_null_model = MCMCglmm(as.formula(formula), 
                                  data = as.data.frame(necropsy_phen_xf), 
                                  random = ~us(trait):Dam + us(trait):NURSE,
                                  rcov = ~us(trait):ID,
                                  prior = prior,
                                  nitt = 103000, burnin = 3000, thin = 100,
                                  family = rep("gaussian", length(necropsy_traits)))

save(necropsy_MCMC_null_model, file = paste0(Rdatas_folder, "necropsy_MCMC_DamNurse_model.Rdata"))
load(paste0(Rdatas_folder, "necropsy_MCMC_DamNurse_model.Rdata"))
load(paste0(Rdatas_folder, "necropsy_MCMC_null_model.Rdata"))
summary(necropsy_MCMC_null_model)

id = colnames(necropsy_MCMC_null_model$VCV)
Gs_dam = array(necropsy_MCMC_null_model$VCV[,grep("Dam", id)], 
                    dim = c(1000, num_necropsy_traits, num_necropsy_traits))
G_dam = apply(array(necropsy_MCMC_null_model$VCV[,grep("Dam", id)], 
                     dim = c(1000, num_necropsy_traits, num_necropsy_traits)), 2:3, median)
G_nurse = apply(array(necropsy_MCMC_null_model$VCV[,grep("NURSE", id)], 
                     dim = c(1000, num_necropsy_traits, num_necropsy_traits)), 2:3, median)
#G_NbyD = apply(array(necropsy_MCMC_null_model$VCV[,grep("NbyD", id)], 
#                      dim = c(1000, num_necropsy_traits, num_necropsy_traits)), 2:3, median)
R_mcmc = apply(array(necropsy_MCMC_null_model$VCV[,grep("ID", id)], 
                     dim = c(1000, num_necropsy_traits, num_necropsy_traits)), 2:3, median)

png("./data/necropsy_fullSib_Dam_Nurse_Gcorrelation.png", width = 2100, height = 700)
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

necropsy_data_cf = filter(necropsy_data, Dam != NURSE)
necropsy_data_ncf = filter(necropsy_data, Dam == NURSE)

stan_cf = stan(file = "./mixedModelGmatrix.stan",
               data = getStanInputMM(necropsy_data_cf, necropsy_traits),
               chain=3, iter = 2000, warmup = 500, control = list(max_treedepth = 11, adapt_delta = 0.99))
stan_ncf = stan(file = "./mixedModelGmatrix.stan",
                data = getStanInputMM(necropsy_data_ncf, necropsy_traits),
                chain=3, iter = 2000, control = list(max_treedepth = 11, adapt_delta = 0.99))

G_cf  = colMeans(rstan::extract(stan_cf, pars = c("G"))[[1]]) * outer(necropsy_phen_sd, necropsy_phen_sd)
G_ncf = colMeans(rstan::extract(stan_ncf, pars = c("G"))[[1]]) * outer(necropsy_phen_sd, necropsy_phen_sd)

png("./data/necropsy_fullSib_CF_nCF_Gcorrelation.png", width = 2100, height = 700)
par(mfrow = c(1, 3), cex=2, oma = c(0, 0, 0, 0))
corrplot.mixed(cov2cor(G_cf),   upper = "ellipse", mar = c(0, 0, 1, 0), main = "FullSib CF")
corrplot.mixed(cov2cor(G_stan), upper = "ellipse", mar = c(0, 0, 1, 0), main = "FullSib Family")
corrplot.mixed(cov2cor(G_ncf),  upper = "ellipse", mar = c(0, 0, 1, 0), main = "FullSib Non-CF")
par(old.par)
dev.off()

save(P_stan, G_stan, Gs_stan, R_stan, G_cf, G_ncf, file = "./Rdatas/necropsy_CovMatrices.Rdata")
load(file = "./Rdatas/necropsy_CovMatrices.Rdata")

### Dam nurse separation
# 
# necropsy_data = inner_join(necropsy_phen, necropsy_markers, by = "ID")
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
#                data = getStanInputDamNurseMM(necropsy_data, necropsy_traits),
#                chain=8, iter = 2000, control = list(max_treedepth = 11, adapt_delta = 0.99))
# save(stan_DamNurse, file = "./Rdatas/necropsy_DamNurseGMatrixFit.Rdata")
# 
# 
# P_DM = colMeans(rstan::extract(stan_DamNurse, pars = c("P"))[[1]])
# G_dam = colMeans(rstan::extract(stan_DamNurse, pars = c("G_dam"))[[1]])
# G_nurse = colMeans(rstan::extract(stan_DamNurse, pars = c("G_nurse"))[[1]])
# Gs_dam = rstan::extract(stan_DamNurse, pars = c("G_dam"))[[1]]
# Gs_nurse = rstan::extract(stan_DamNurse, pars = c("G_nurse"))[[1]]
# R_DM = colMeans(rstan::extract(stan_DamNurse, pars = c("R"))[[1]])
# save(P_DM, G_dam, G_nurse, Gs_dam, Gs_nurse, R_DM, file = "./Rdatas/necropsy_DamNurseMatrices.Rdata")
# load(file = "./Rdatas/necropsy_DamNurseMatrices.Rdata")
# 
# 
# par(mfrow = c(1, 3), cex=2, oma = c(0, 0, 0, 0))
# corrplot.mixed(cov2cor(G_dam),   upper = "ellipse", mar = c(0, 0, 1, 0), main = "FullSib Dam")
# corrplot.mixed(cov2cor(G_stan), upper = "ellipse", mar = c(0, 0, 1, 0), main = "FullSib Family")
# corrplot.mixed(cov2cor(G_nurse), upper = "ellipse", mar = c(0, 0, 1, 0), main = "FullSib Nurse")



