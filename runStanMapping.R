setwd("/home/diogro/projects/mouse-qtls")
source('read_mouse_data.R')
source('utils.R')

#Rdatas_folder = "~/gdrive/LGSM_project_Rdatas/"
Rdatas_folder = "./data/Rdatas/"

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = 80)

area_data = inner_join(area_phen_std,
                       Reduce(inner_join, area_markers),
                       by = "ID")

area_data = area_data %>% mutate_each(funs(scale), matches('area'))

ddply(area_data, .(FAMILY), function(x) select(x, matches('area')) %>% colMeans) %>% summary

current_chrom = 6
getStanInput = function(current_chrom){
    K        = num_area_traits
    J        = loci_per_chrom[current_chrom]
    N        = dim(area_data)[1]
    n_family = length(unique(area_data$FAMILY))
    family   = as.integer(as.factor(area_data$FAMILY))
    ad       = as.matrix(select(area_data, matches('chrom6_A')))
    dm       = as.matrix(select(area_data, matches('chrom6_D')))
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
stan_parameters = getStanInput(6)
names(stan_parameters)

stan_model_SUR_HC = stan(file = './SUR_horseShoe.stan',
                         data = stan_parameters, chain=2, iter = 100)

getStanEffects = function(stan_model){
  HC_summary = summary(stan_model, pairs = c("w_ad", "w_dm"))$summary
  s = loci_per_chrom[current_chrom] * num_area_traits * 2
  mask = grepl("w_", rownames(HC_summary))
  effects = data.frame(HC_summary[mask, c("50%", "2.5%", "97.5%")]) 
  colnames(effects) <- c("median", "lower", "upper")
  effects$type = rep(c("additive", "dominance"), each = s/2)
  effects$chrom = current_chrom
  effects$marker = rep(1:loci_per_chrom[current_chrom], 2*num_area_traits)
  effects$trait = rep(area_traits, each = loci_per_chrom[current_chrom])
  tbl_df(effects)
}

effects = getStanEffects(stan_model_SUR_HC)
hc_plot = ggplot(effects, aes(marker, median, group = trait)) +
  geom_point() + facet_grid(trait~type) +
  geom_hline(yintercept = 0) +
  geom_point(size = 0.3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, size = 0.3)
save_plot("data/figures/stan_SUR_HC_chrom6.png", hc_plot, base_height = 6, base_aspect_ratio = 1.8)
