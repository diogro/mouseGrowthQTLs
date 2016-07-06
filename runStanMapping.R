setwd("/home/diogro/projects/mouse-qtls")
source('read_mouse_data.R')
source('utils.R')

#Rdatas_folder = "~/gdrive/LGSM_project_Rdatas/"
Rdatas_folder = "./data/Rdatas/"

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = 80)

area_data = inner_join(area_phen_std,
                       simulated_genomes[[4]],
                       by = "ID")

area_data = area_data %>% mutate_each(funs(scale), matches('area'))

ddply(area_data, .(FAMILY), function(x) select(x, matches('area')) %>% colMeans) %>% summary

true_effects = c(findA(area_data$area1, area_data$chrom4_A10, 0.1),
                 findA(area_data$area2, area_data$chrom4_A10, 0.1),
                 findA(area_data$area3, area_data$chrom4_A13, 0.1),
                 findA(area_data$area4, area_data$chrom4_A14, 0.1),
                 findA(area_data$area5, area_data$chrom4_A14, 0.1),
                 findA(area_data$area6, area_data$chrom4_A13, 0.1),
                 findA(area_data$area7, area_data$chrom4_A14, 0.1))
true_effects = data.frame(true_effects, trait = area_traits, type = "additive",
                          marker = c(rep(10, 2), 13, rep(14, 2), 13, 14))

area_data$area1 = makeSimData(area_data$area1, area_data$chrom4_A10, 0.1)
area_data$area2 = makeSimData(area_data$area2, area_data$chrom4_A10, 0.1)
area_data$area3 = makeSimData(area_data$area3, area_data$chrom4_A13, 0.1)
area_data$area4 = makeSimData(area_data$area4, area_data$chrom4_A14, 0.1)
area_data$area5 = makeSimData(area_data$area5, area_data$chrom4_A14, 0.1)
area_data$area6 = makeSimData(area_data$area6, area_data$chrom4_A13, 0.1)
area_data$area7 = makeSimData(area_data$area7, area_data$chrom4_A14, 0.1)

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
stan_model_G = stan(file = './mixedModelGmatrix.stan',
                    data = stan_parameters[-c(6,7)], chain=4, iter = 1000)

stan_model_SUR_HC = stan(file = './SUR_horseShoe.stan',
                         data = stan_parameters, chain=4, iter = 200)

HC_model = stan_model(file = './SUR_horseShoe.stan')
vb_HC_model = vb(HC_model, data = stan_parameters)

stan_model = stan_model_SUR_HC
stan_model = vb_HC_model
getStanEffects = function(stan_model){
  HC_summary = summary(stan_model, pairs = c("w_ad", "w_dm"))$summary
  mask = grepl("w_", rownames(HC_summary))
  effects = data.frame(HC_summary[mask, c("mean", "2.5%", "97.5%")]) 
  colnames(effects) <- c("mean", "lower", "upper")
  effects$type = rep(c("additive", "dominance"), each = s/2)
  effects$chrom = current_chrom
  effects$marker = rep(1:loci_per_chrom[current_chrom], 2*num_area_traits)
  effects$trait = rep(area_traits, each = loci_per_chrom[current_chrom])
  tbl_df(effects)
}

effects = getStanEffects(vb_HC_model)
hc_plot = ggplot(effects, aes(marker, mean, group = trait)) +
  geom_point() + facet_grid(trait~type, scales = "free") +
  geom_hline(yintercept = 0) +
  geom_point(size = 0.3) +
  geom_point(data = true_effects, aes(y = true_effects), color = "red") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, size = 0.3)
save_plot("data/figures/sim_stan_SUR_HC_vb.png", hc_plot, base_height = 6, base_aspect_ratio = 1.8)

effects = getStanEffects(stan_model_SUR_HC)
hc_plot = ggplot(effects, aes(marker, mean, group = trait)) +
  geom_point() + facet_grid(trait~type, scales = "free") +
  geom_hline(yintercept = 0) +
  geom_point(size = 0.3) +
  geom_point(data = true_effects, aes(y = true_effects), color = "red") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, size = 0.3)
save_plot("data/figures/sim_stan_SUR_HC.png", hc_plot, base_height = 6, base_aspect_ratio = 1.8)
