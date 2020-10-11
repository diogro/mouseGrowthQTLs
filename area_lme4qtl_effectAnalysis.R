setwd("/home/diogro/projects/mouse-qtls")
source('read_mouse_data.R')

#Rdatas_folder = "~/gdrive/LGSM_project_Rdatas/"
Rdatas_folder = "./Rdatas/"

if(!require(doMC)){install.packages("doMC"); library(doMC)}
if(!require(lme4qtl)){install.packages("lme4qtl"); library(lme4qtl)}
if(!require(qvalue)){install.packages("qvalue"); library(qvalue)}
registerDoMC(3)

area_data = inner_join(area_phen_std,
                       area_markers,
                       by = "ID") %>%
  gather(variable, value, area1:area7)
area_data = mutate(area_data, FAMILY = as.factor(FAMILY), variable = as.factor(variable))

e_SUR = read_csv("./data/area_significant_marker_effects_SUR.csv") %>% arrange(trait, id, class)
e_lm = read_csv("./data/area_significant_marker_effects_lm.csv") %>% arrange(trait, id, class)
load(file = "./Rdatas/area_significant_stan_fit.Rdata")
e_HC = read_csv("./data/area_significant_marker_effects_HC.csv") %>% arrange(trait, id, class)
load(file = "./Rdatas/area_significant_stan_HCp_fit.Rdata")
full_HCp = readRDS("./Rdatas/area_scaled_allmarkers_HCPlus")

effect_matrix_lm = e_lm %>%
  select(id, class, trait, mean) %>%
  spread(trait, mean) %>%
  filter(class == "additive") %>% select(-class)
effect_matrix_SUR = e_SUR %>%
  select(id, class, trait, mean) %>%
  spread(trait, mean) %>%
  filter(class == "additive") %>% select(-class)
effect_matrix_HC = e_HC %>%
  select(id, class, trait, mean) %>%
  spread(trait, mean) %>%
  filter(class == "additive") %>% select(-class)
png("data/area_additive_effects_PCA.png")
biplot(prcomp(effect_matrix_SUR[,area_traits]))
dev.off()
eig_effects  =  eigen(cov(effect_matrix_SUR[,area_traits]))
plot(eig_effects$values)
cov2cor(G)
plot(e_SUR$mean, e_HC$mean)
abline(0, 1)

source("./area_parental_strains.R")
LG = filter(ancestral_areas, strain == "LG") %>% select(area_traits) %>% as.numeric
SM = filter(ancestral_areas, strain == "SM") %>% select(area_traits) %>% as.numeric
F3 = sapply(area_phen[,area_traits], mean)

d_z = LG - SM
dz = LG - SM
library(evolqg)
area_sds = apply(area_phen[,area_traits], 2, sd)
G = G_stan * outer(area_sds, area_sds)
corrplot.mixed(cov2cor(G), upper = "ellipse")
plot(eigen(G_mcmc)$values)
G_ext4 = ExtendMatrix(G, ret.dim = 4)[[1]]
G_ext5 = ExtendMatrix(G, ret.dim = 5)[[1]]
solve(G, d_z)
solve(G_ext4, d_z)
beta = solve(G_ext5, d_z)

library(viridis)
vectorCor = function(x, y) Normalize(x) %*% Normalize(y)
vectorCor(d_z, beta)
random_vec = matrix(rnorm(7*1000), 1000, 7)
quantile(abs(apply(random_vec, 1, vectorCor, rep(1, 7))), 0.95)
crss = data.frame(beta = apply(effect_matrix[,growth_traits], 1, vectorCor, beta),
                  dz = apply(effect_matrix[,growth_traits], 1, vectorCor, d_z)) %>% gather(class, value, beta:dz)

corrs = data.frame(beta = apply(effect_matrix[,growth_traits], 1, vectorCor, beta),
                  dz = apply(effect_matrix[,growth_traits], 1, vectorCor, d_z),
                  norm = apply(effect_matrix[,growth_traits], 1, Norm))
ggplot(crss, aes(class, value, fill = class)) + geom_violin()
ggplot(corrs, aes(norm, beta)) + geom_point() + geom_smooth(method = "lm")
ggplot(corrs, aes(norm, dz)) + geom_point() + geom_smooth(method = "lm")

area_m = as.numeric(ddply(area_phen, .(SEX), numcolwise(mean))[2,area_traits])
area_f = as.numeric(ddply(area_phen, .(SEX), numcolwise(mean))[1,area_traits])

area_prediction = reshape2::melt(data.frame(trait = as.factor(area_traits), 
                                            SM_Observed = SM,
                                            F3_Observed = F3,
                                            LG_Observed = LG)) %>% separate(variable, c("Line", "Type"))
area_pred_plot = ggplot(area_prediction, aes(trait, value, group = interaction(Type, Line), color = Line, linetype = Type)) + geom_line(size = 1) + scale_x_discrete(labels = paste("area", 1:7)) + labs(y = "Area (mm2)", x = "Start area")


posterior_predict = function(beta_ad){
  SM_e = rowSums(-1 * beta_ad) * sapply(area_phen[,area_traits], sd) + 
    sapply(area_phen[,area_traits], mean)
  LG_e = rowSums(beta_ad) * sapply(area_phen[,area_traits], sd) + 
    sapply(area_phen[,area_traits], mean)
  reshape2::melt(data.frame(trait = as.factor(area_traits), 
                            SM_Predicted = SM_e,
                            LG_Predicted = LG_e)) %>% separate(variable, c("Line", "Type"))
}

w_ad = rstan::extract(stan_model_SUR, pars = "w_ad")[[1]]
effect_matrix = aaply(w_ad, c(2, 3), mean)
SM_e = rowSums(-1 * as.matrix(effect_matrix)) * sapply(area_phen[,area_traits], sd) + 
  sapply(area_phen[,area_traits], mean)
LG_e = rowSums(as.matrix(effect_matrix)) * sapply(area_phen[,area_traits], sd) + 
  sapply(area_phen[,area_traits], mean)
post = adply(w_ad, 1, posterior_predict)
area_prediction = reshape2::melt(data.frame(trait = as.factor(area_traits), 
                                            SM_stan = SM_e,
                                            SM_Observed = SM,
                                            LG_stan = LG_e,
                                            F3_Observed = F3,
                                            LG_Observed = LG)) %>% separate(variable, c("Line", "Type"))
area_pred_plot_SUR = ggplot() + scale_x_discrete(labels = paste("area", 1:7)) + labs(y = "Area (mm2)", x = "Start area") + geom_line(data = post, color = "gray", size = 0.5, linetype = 1, alpha = 0.1, aes(trait, value, group = interaction(Type, Line, iterations))) + geom_line(size = 1, data = area_prediction, aes(trait, value, group = interaction(Type, Line), color = Line, linetype = Type))
save_plot("data/area_multiple_regression_ancestral_prediction.png", area_pred_plot_SUR, base_height = 7, base_aspect_ratio = 2)

w_ad = rstan::extract(full_HCp, pars = "w_ad")[[1]]
effect_matrix = aaply(w_ad, c(2, 3), mean)
SM_e = rowSums(-1 * as.matrix(effect_matrix)) * sapply(area_phen[,area_traits], sd) + 
  sapply(area_phen[,area_traits], mean)
LG_e = rowSums(as.matrix(effect_matrix)) * sapply(area_phen[,area_traits], sd) + 
  sapply(area_phen[,area_traits], mean)
post = adply(w_ad, 1, posterior_predict)
area_prediction = reshape2::melt(data.frame(trait = as.factor(area_traits), 
                                            SM_stan = SM_e,
                                            SM_Observed = SM,
                                            LG_stan = LG_e,
                                            F3_Observed = F3,
                                            LG_Observed = LG)) %>% separate(variable, c("Line", "Type"))
area_pred_plot_HC_full = ggplot() + scale_x_discrete(labels = paste("area", 1:7)) + labs(y = "Area (mm2)", x = "Start area") + geom_line(data = post, color = "gray", size = 0.5, linetype = 1, alpha = 0.1, aes(trait, value, group = interaction(Type, Line, iterations))) + geom_line(size = 1, data = area_prediction, aes(trait, value, group = interaction(Type, Line), color = Line, linetype = Type))
plot_grid(area_pred_plot_SUR, area_pred_plot_HC_full)

effect_matrix[,1]
Reduce("+", alply(effect_matrix, 2, function(x) outer(x, x)))
