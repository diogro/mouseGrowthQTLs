setwd("/home/diogro/projects/mouse-qtls")
source('read_mouse_data.R')
source('utils.R')

Rdatas_folder = "~/gdrive/LGSM_project_Rdatas/"
#Rdatas_folder = "./data/Rdatas/"

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = 3)

area_data = inner_join(area_phen_std,
                       simulated_genomes[[5]],
                       by = "ID")

area_data[area_traits] = scale(area_data[area_traits])
area_data[area_traits] = scale(area_data[area_traits], scale = rep(sqrt(1/20), num_area_traits))

ddply(area_data, .(FAMILY), function(x) select(x, matches('area')) %>% colMeans)

true_effects = c(findA(area_data$area1, area_data$chrom4_A10, 0.01),
                 findA(area_data$area2, area_data$chrom4_A10, 0.01),
                 findA(area_data$area3, area_data$chrom4_A13, 0.01),
                 findA(area_data$area4, area_data$chrom4_A14, 0.01),
                 findA(area_data$area5, area_data$chrom4_A14, 0.01),
                 findA(area_data$area6, area_data$chrom4_A13, 0.01),
                 findA(area_data$area7, area_data$chrom4_A14, 0.01))
true_effects = data.frame(true_effects, trait = area_traits, type = "additive",
                          marker = c(rep(10, 2), 13, rep(14, 2), 13, 14), chrom = 4)

area_data$area1 = makeSimData(area_data$area1, area_data$chrom4_A10, 0.01)
area_data$area2 = makeSimData(area_data$area2, area_data$chrom4_A10, 0.01)
area_data$area3 = makeSimData(area_data$area3, area_data$chrom4_A13, 0.01)
area_data$area4 = makeSimData(area_data$area4, area_data$chrom4_A14, 0.01)
area_data$area5 = makeSimData(area_data$area5, area_data$chrom4_A14, 0.01)
area_data$area6 = makeSimData(area_data$area6, area_data$chrom4_A13, 0.01)
area_data$area7 = makeSimData(area_data$area7, area_data$chrom4_A14, 0.01)

#area_data = area_data %>% mutate_each(funs(scale), matches('area'))
#area_data[area_traits] = scale(area_data[area_traits], scale = rep(0.1, num_area_traits))

ddply(area_data, .(FAMILY), function(x) select(x, matches('area')) %>% colMeans) %>% summary

current_chrom = 4
getStanInput = function(current_chrom){
    K        = num_area_traits
    J        = loci_per_chrom[current_chrom]
    N        = dim(area_data)[1]
    n_family = length(unique(area_data$FAMILY))
    family   = as.integer(as.factor(area_data$FAMILY))
    ad       = as.matrix(select(area_data, matches('chrom4_A')))
    dm       = as.matrix(select(area_data, matches('chrom4_D')))
    y        = as.matrix(select(area_data, matches('area')))
    beta_ad  = matrix(0., K, J)
    beta_dm  = matrix(0., K, J)
    param_list = list(K        = K,
                      J        = J,
                      N        = N,
                      n_family = n_family,
                      family   = family,
                      ad       = ad,
                      dm       = dm,
                      y        = y,
                      beta_ad  = beta_ad,
                      beta_dm  = beta_dm)
    return(param_list)
}
stan_parameters = getStanInput(4)
names(stan_parameters)

stan_model_SUR_HC <- 0
stan_model_SUR_HC = stan(file = './SUR_horseShoe.stan',
                         data = stan_parameters, chain=3, iter = 300)

x = diag(summary(stan_model_SUR_HC, pars = "L_sigma_G")[[1]][,1]) %*%
 matrix(summary(stan_model_SUR_HC, pars = "L_Omega_G")[[1]][,1], 7, 7, T)
G = x %*% t(x)
x = diag(summary(stan_model_SUR_HC, pars = "L_sigma_R")[[1]][,1]) %*%
 matrix(summary(stan_model_SUR_HC, pars = "L_Omega_R")[[1]][,1], 7, 7, T)
R = x %*% t(x)
P = area_data %>% select(matches('area')) %>% cov
g_p = data.frame(P = P[lower.tri(P, diag = T)], G_R = (G+R)[lower.tri(P, diag = T)])
gp_plot = ggplot(g_p, aes(P, G_R)) + geom_point() + geom_abline()
#save_plot("./data/figures/gp_plot.png", gp_plot, base_height = 6, base_aspect_ratio = 1.8)

shrink = (extract(stan_model_SUR_HC, pars = c("lambda_ad", "tau")))
weigths = sweep(shrink[[1]], 1, shrink[[2]], "*")
reshape2::melt(weigths) %>% ddply(.(Var2, Var3), numcolwise(mean)) %>%
    ggplot(aes(Var3, value)) + geom_point() + geom_hline(yintercept = 0.2) + facet_wrap(~Var2, scale = "free_x")
plot(stan_model_SUR_HC, pars = c("beta_ad"))
plot(stan_model_SUR_HC, pars = c("lambda_ad"))
plot(stan_model_SUR_HC, pars = c("w_ad"))
plot(stan_model_SUR_HC, pars = c("tau"))

stan_model_SUR_HCPlus = stan(file = './SUR_horseShoePlus.stan',
                         data = stan_parameters, chain=10, iter = 200)

getStanEffects = function(stan_model){
  HC_summary = summary(stan_model, pairs = c("w_ad", "w_dm"))$summary
  s = loci_per_chrom[current_chrom] * num_area_traits * 2
  mask = grepl("w_", rownames(HC_summary))
  effects = data.frame(HC_summary[mask, c("mean", "2.5%", "97.5%")])
  colnames(effects) <- c("mean", "lower", "upper")
  effects$type = rep(c("aditive", "dominance"), each = s/2)
  effects$chrom = current_chrom
  effects$marker = rep(1:loci_per_chrom[current_chrom], 2*num_area_traits)
  effects$trait = rep(area_traits, each = loci_per_chrom[current_chrom])
  tbl_df(effects)
}
getStanShrinkage = function(stan_model){
  shrink = (extract(stan_model, pars = c("lambda_ad", "lambda_dm", "tau")))
  raw_weights_ad = sweep(shrink[[1]], 1, shrink[[3]], "*")
  raw_weights_dm = sweep(shrink[[2]], 1, shrink[[3]], "*")
  weights_ad = reshape2::melt(raw_weights_ad) %>% ddply(.(Var2, Var3), numcolwise(mean))
  weights_dm = reshape2::melt(raw_weights_dm) %>% ddply(.(Var2, Var3), numcolwise(mean))
  weights = rbind(select(weights_ad, -iterations), select(weights_dm, -iterations))
  s = loci_per_chrom[current_chrom] * num_area_traits * 2
  colnames(weights) <- c("trait", "marker", "mean")
  weights$type = rep(c("aditive", "dominance"), each = s/2)
  weights$chrom = current_chrom
  weights$marker = rep(1:loci_per_chrom[current_chrom], 2*num_area_traits)
  weights$trait = rep(area_traits, each = loci_per_chrom[current_chrom])
  tbl_df(weights)
}

plotEffectEstimate = function(effects, current_chrom, file = NULL){
    hc_plot = ggplot(filter(effects, chrom == current_chrom), aes(marker, mean, group = trait)) +
        geom_point() + facet_grid(trait~type, scales = "free") +
        geom_hline(yintercept = 0) +
        geom_point(size = 0.3) +
        geom_point(data = filter(true_effects, chrom == current_chrom), aes(y = true_effects), color = "red") +
        geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, size = 0.3)
    if(!is.null(file))
        save_plot(file, hc_plot, base_height = 6, base_aspect_ratio = 1.8)
    return(hc_plot)
}
plotWeights = function(weights, current_chrom, file = NULL){
    hc_plot = ggplot(filter(weights, chrom == current_chrom), aes(marker, mean, group = trait)) +
        geom_point() + facet_grid(trait~type, scales = "free") +
        geom_hline(yintercept = 0) +
        geom_point(size = 0.3)
    if(!is.null(file))
        save_plot(file, hc_plot, base_height = 6, base_aspect_ratio = 1.8)
    return(hc_plot)
}

effects = getStanEffects(stan_model_SUR_HC)
plotEffectEstimate(effects, 4, file=  "./data/figures/sim_chrom4_StanSUR_2pct.png")
weights = getStanShrinkage(stan_model_SUR_HC)
plotWeights(weights, 4)

effects = getStanEffects(stan_model_SUR_HCPlus)
plotEffectEstimate(4, file=  "./data/figures/sim_chrom4_StanSURPlus_3pct.png")
