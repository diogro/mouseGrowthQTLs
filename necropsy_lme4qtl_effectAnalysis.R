setwd("/home/diogro/projects/mouse-qtls")
source('read_mouse_data.R')

#Rdatas_folder = "~/gdrive/LGSM_project_Rdatas/"
Rdatas_folder = "./Rdatas/"

install_load("doMC", "lme4qtl", "qvalue")
registerDoMC(2)

necropsy_data = inner_join(necropsy_phen_std,
                           necropsy_markers,
                           by = "ID") %>%
gather(variable, value, FATPAD:SPLEEN)
necropsy_data = mutate(necropsy_data, FAMILY = as.factor(FAMILY), variable = as.factor(variable))

e_SUR = read_csv("./data/necropsy_significant_marker_effects_SUR.csv") %>% arrange(trait, id, class)
e_lm = read_csv("./data/necropsy_significant_marker_effects_lm.csv") %>% arrange(trait, id, class) %>% filter(class != "imprinting")
load(file = "./Rdatas/necropsy_significant_stan_fit.Rdata")
e_HC = read_csv("./data/necropsy_significant_marker_effects_SUR_HC.csv") %>% arrange(trait, id, class)
load(file = "./Rdatas/necropsy_significant_stan_HCp_fit.Rdata")
full_HCp = readRDS("./Rdatas/necropsy_scaled_allmarkers_HCPlus")

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
png("data/necropsy_additive_effects_PCA.png")
par(mfrow = c(1, 1))
biplot(prcomp(effect_matrix_SUR[,necropsy_traits]))
dev.off()
eig_effects  =  eigen(cov(effect_matrix_lm[,necropsy_traits]))
plot(eig_effects$values)
cov2cor(G)
plot(e_SUR$mean, e_lm$mean)
abline(0, 1)

source("./necropsy_parental_strains.R")
LG = filter(ancestral_necropsys, strain == "LG") %>% select(necropsy_traits) %>% as.numeric
SM = filter(ancestral_necropsys, strain == "SM") %>% select(necropsy_traits) %>% as.numeric
F3 = sapply(necropsy_phen[,necropsy_traits], mean)


necropsy_m = as.numeric(ddply(necropsy_phen, .(SEX), numcolwise(mean))[2,necropsy_traits])
necropsy_f = as.numeric(ddply(necropsy_phen, .(SEX), numcolwise(mean))[1,necropsy_traits])

necropsy_prediction = reshape2::melt(data.frame(trait = as.factor(necropsy_traits), 
                                                SM_Observed = SM,
                                                F3_Observed = F3,
                                                LG_Observed = LG)) %>% separate(variable, c("Line", "Type"))
necropsy_pred_plot = ggplot(necropsy_prediction, aes(trait, value, group = interaction(Type, Line), color = Line, linetype = Type)) + geom_line(size = 1) + scale_x_discrete(labels = paste("necropsy", 1:7)) + labs(y = "necropsy (mm2)", x = "Start necropsy")


posterior_predict = function(beta_ad){
    SM_e = rowSums(-1 * beta_ad) * sapply(necropsy_phen[,necropsy_traits], sd) + 
        sapply(necropsy_phen[,necropsy_traits], mean)
    LG_e = rowSums(beta_ad) * sapply(necropsy_phen[,necropsy_traits], sd) + 
        sapply(necropsy_phen[,necropsy_traits], mean)
    reshape2::melt(data.frame(trait = as.factor(necropsy_traits), 
                              SM_Predicted = SM_e,
                              LG_Predicted = LG_e)) %>% separate(variable, c("Line", "Type"))
}

w_ad = rstan::extract(stan_model_SUR, pars = "w_ad")[[1]]
effect_matrix = aaply(w_ad, c(2, 3), mean)
SM_e = rowSums(-1 * as.matrix(effect_matrix)) * sapply(necropsy_phen[,necropsy_traits], sd) + 
    sapply(necropsy_phen[,necropsy_traits], mean)
LG_e = rowSums(as.matrix(effect_matrix)) * sapply(necropsy_phen[,necropsy_traits], sd) + 
    sapply(necropsy_phen[,necropsy_traits], mean)
post = adply(w_ad, 1, posterior_predict)
necropsy_prediction = reshape2::melt(data.frame(trait = as.factor(necropsy_traits), 
                                                SM_stan = SM_e,
                                                SM_Observed = SM,
                                                LG_stan = LG_e,
                                                F3_Observed = F3,
                                                LG_Observed = LG)) %>% separate(variable, c("Line", "Type"))
necropsy_pred_plot_SUR = ggplot() + scale_x_discrete(labels = paste("necropsy", 1:7)) + labs(y = "necropsy (mm2)", x = "Start necropsy") + geom_line(data = post, color = "gray", size = 0.5, linetype = 1, alpha = 0.1, aes(trait, value, group = interaction(Type, Line, iterations))) + geom_line(size = 1, data = necropsy_prediction, aes(trait, value, group = interaction(Type, Line), color = Line, linetype = Type))
save_plot("data/necropsy_multiple_regression_ancestral_prediction.png", necropsy_pred_plot_SUR, base_height = 7, base_aspect_ratio = 2)

w_ad = rstan::extract(full_HCp, pars = "w_ad")[[1]]
effect_matrix = aaply(w_ad, c(2, 3), mean)
SM_e = rowSums(-1 * as.matrix(effect_matrix)) * sapply(necropsy_phen[,necropsy_traits], sd) + 
    sapply(necropsy_phen[,necropsy_traits], mean)
LG_e = rowSums(as.matrix(effect_matrix)) * sapply(necropsy_phen[,necropsy_traits], sd) + 
    sapply(necropsy_phen[,necropsy_traits], mean)
post = adply(w_ad, 1, posterior_predict)
necropsy_prediction = reshape2::melt(data.frame(trait = as.factor(necropsy_traits), 
                                                SM_stan = SM_e,
                                                SM_Observed = SM,
                                                LG_stan = LG_e,
                                                F3_Observed = F3,
                                                LG_Observed = LG)) %>% separate(variable, c("Line", "Type"))
necropsy_pred_plot_HC_full = ggplot() + scale_x_discrete(labels = paste("necropsy", 1:7)) + labs(y = "necropsy (mm2)", x = "Start necropsy") + geom_line(data = post, color = "gray", size = 0.5, linetype = 1, alpha = 0.1, aes(trait, value, group = interaction(Type, Line, iterations))) + geom_line(size = 1, data = necropsy_prediction, aes(trait, value, group = interaction(Type, Line), color = Line, linetype = Type))
plot_grid(necropsy_pred_plot_SUR, necropsy_pred_plot_HC_full)

effect_matrix[,1]
Reduce("+", alply(effect_matrix, 2, function(x) outer(x, x)))
