setwd("/home/diogro/projects/mouse-qtls")
source('read_mouse_data.R')

#Rdatas_folder = "~/gdrive/LGSM_project_Rdatas/"
Rdatas_folder = "./Rdatas/"

install_load("doMC", "lme4qtl", "qvalue")
registerDoMC(10)

growth_data = inner_join(growth_phen_std,
                        growth_markers,
                        by = "ID") %>%
  gather(variable, value, growth12:growth78)
growth_data = mutate(growth_data, FAMILY = as.factor(FAMILY), variable = as.factor(variable))

e = read_csv("./data/growth_significant_marker_effects_SUR.csv") %>% arrange(trait, id, class)
load(file = "./Rdatas/significant_stan_fit.Rdata")
eHC = read_csv("./data/growth_significant_marker_effectsHC.csv") %>% arrange(trait, id, class)
load(file = "./Rdatas/significant_stan_HCp_fit.Rdata")
full_HCp = readRDS("./Rdatas/growth_scaled_allmarkers_HCPlus")

effect_matrix = e %>%
    select(id, class, trait, mean) %>%
    spread(trait, mean) %>%
    filter(class == "dominance") %>% select(-class)
effect_matrix_HC = eHC %>%
    select(id, class, trait, mean) %>%
    spread(trait, mean) %>%
    filter(class == "additive") %>% select(-class)
png("data/growth_additive_effects_PCA.png")
biplot(prcomp(effect_matrix[,growth_traits]))
dev.off()
eig_effects  =  eigen(cov(effect_matrix[,growth_traits]))
plot(eig_effects$values)
cov2cor(G)
plot(e$mean, eHC$mean)
abline(0, 1)


LG = c(3.785,4.435,8.43,7.395,2.995,1.85,2.085)
SM = c(3.31 ,2.98,3.82,2.175,0.765,1.165,0.51)
F3 = sapply(growth_phen[,growth_traits], mean)

d_z = LG - SM
dz = LG - SM
library(evolqg)
growth_sds = apply(growth_phen[,growth_traits], 2, sd)
G = G_stan * outer(growth_sds, growth_sds)
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



growth_m = as.numeric(ddply(growth_phen, .(SEX), numcolwise(mean))[2,growth_traits])
growth_f = as.numeric(ddply(growth_phen, .(SEX), numcolwise(mean))[1,growth_traits])

growth_prediction = reshape2::melt(data.frame(trait = as.factor(growth_traits), 
                                              SM_SUR = SM_e,
                                              SM_HC = SM_e_HC,
                                              SM_Observed = SM,
                                              LG_SUR = LG_e,
                                              LG_HC = LG_e_HC,
                                              F3_Observed = F3,
                                              LG_Observed = LG)) %>% separate(variable, c("Line", "Type"))
growth_pred_plot = ggplot(growth_prediction, aes(trait, value, group = interaction(Type, Line), color = Line, linetype = Type)) + geom_line(size = 1) + scale_x_discrete(labels = paste("Week", 1:7)) + labs(y = "Weekly growth (g)", x = "Start week")


posterior_predict = function(beta_ad){
    SM_e = rowSums(-1 * beta_ad) * sapply(growth_phen[,growth_traits], sd) + 
        sapply(growth_phen[,growth_traits], mean)
    LG_e = rowSums(beta_ad) * sapply(growth_phen[,growth_traits], sd) + 
        sapply(growth_phen[,growth_traits], mean)
    reshape2::melt(data.frame(trait = as.factor(growth_traits), 
                                              SM_Predicted = SM_e,
                                              LG_Predicted = LG_e)) %>% separate(variable, c("Line", "Type"))
}

w_ad = rstan::extract(stan_model_SUR, pars = "w_ad")[[1]]
effect_matrix = aaply(w_ad, c(2, 3), mean)
SM_e = rowSums(-1 * as.matrix(effect_matrix)) * sapply(growth_phen[,growth_traits], sd) + 
    sapply(growth_phen[,growth_traits], mean)
LG_e = rowSums(as.matrix(effect_matrix)) * sapply(growth_phen[,growth_traits], sd) + 
    sapply(growth_phen[,growth_traits], mean)
post = adply(w_ad, 1, posterior_predict)
growth_prediction = reshape2::melt(data.frame(trait = as.factor(growth_traits), 
                                              SM_stan = SM_e,
                                              SM_Observed = SM,
                                              LG_stan = LG_e,
                                              F3_Observed = F3,
                                              LG_Observed = LG)) %>% separate(variable, c("Line", "Type"))
growth_pred_plot_SUR = ggplot() + scale_x_discrete(labels = paste("Week", 1:7)) + labs(y = "Weekly growth (g)", x = "Start week") + geom_line(data = post, color = "gray", size = 0.5, linetype = 1, alpha = 0.1, aes(trait, value, group = interaction(Type, Line, iterations))) + geom_line(size = 1, data = growth_prediction, aes(trait, value, group = interaction(Type, Line), color = Line, linetype = Type))
save_plot("data/growth_multiple_regression_ancestral_prediction.png", growth_pred_plot_SUR, base_height = 7, base_aspect_ratio = 2)

w_ad = rstan::extract(full_HCp, pars = "w_ad")[[1]]
effect_matrix = aaply(w_ad, c(2, 3), mean)
SM_e = rowSums(-1 * as.matrix(effect_matrix)) * sapply(growth_phen[,growth_traits], sd) + 
    sapply(growth_phen[,growth_traits], mean)
LG_e = rowSums(as.matrix(effect_matrix)) * sapply(growth_phen[,growth_traits], sd) + 
    sapply(growth_phen[,growth_traits], mean)
post = adply(w_ad, 1, posterior_predict)
growth_prediction = reshape2::melt(data.frame(trait = as.factor(growth_traits), 
                                              SM_stan = SM_e,
                                              SM_Observed = SM,
                                              LG_stan = LG_e,
                                              F3_Observed = F3,
                                              LG_Observed = LG)) %>% separate(variable, c("Line", "Type"))
growth_pred_plot_HC_full = ggplot() + scale_x_discrete(labels = paste("Week", 1:7)) + labs(y = "Weekly growth (g)", x = "Start week") + geom_line(data = post, color = "gray", size = 0.5, linetype = 1, alpha = 0.1, aes(trait, value, group = interaction(Type, Line, iterations))) + geom_line(size = 1, data = growth_prediction, aes(trait, value, group = interaction(Type, Line), color = Line, linetype = Type))
plot_grid(growth_pred_plot_SUR, growth_pred_plot_HC_full)

effect_matrix[,1]
Reduce("+", alply(effect_matrix, 2, function(x) outer(x, x)))
