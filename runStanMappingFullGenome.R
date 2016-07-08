setwd("/home/diogro/projects/mouse-qtls")
source('read_mouse_data.R')
source('utils.R')

#Rdatas_folder = "~/gdrive/LGSM_project_Rdatas/"
Rdatas_folder = "./data/Rdatas/"

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = 20)

area_data = inner_join(area_phen_std,
                       Reduce(inner_join, area_markers),
                       by = "ID")

area_data = area_data %>% mutate_each(funs(scale), matches('area'))

ddply(area_data, .(FAMILY), function(x) select(x, matches('area')) %>% colMeans)

getStanInput = function(){
    K        = num_area_traits
    J        = sum(loci_per_chrom)
    N        = dim(area_data)[1]
    n_family = length(unique(area_data$FAMILY))
    family   = as.integer(as.factor(area_data$FAMILY))
    ad       = as.matrix(select(area_data, matches('_A')))
    dm       = as.matrix(select(area_data, matches('_D')))
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
stan_parameters = getStanInput()
names(stan_parameters)

stan_model_SUR_HC = stan(file = './SUR_horseShoe.stan',
                         data = stan_parameters, chain=20, iter = 100)

stan_model = stan_model_SUR_HC
getStanEffects = function(stan_model){
  HC_summary = summary(stan_model, pairs = c("w_ad", "w_dm"))$summary
  s = sum(loci_per_chrom) * num_area_traits * 2
  mask = grepl("w_", rownames(HC_summary))
  effects = data.frame(HC_summary[mask, c("mean", "2.5%", "97.5%")])
  colnames(effects) <- c("mean", "lower", "upper")
  effects$type = rep(c("additive", "dominance"), each = s/2)
  effects$chrom = rep(unlist(lapply(seq_along(loci_per_chrom), function(x) rep(x, loci_per_chrom[x]))), 2)
  effects$marker = rep(unlist(lapply(loci_per_chrom, function(x) 1:x)), 2)
  effects$trait = rep(area_traits, each = sum(loci_per_chrom))
  tbl_df(effects)
}

effects = getStanEffects(stan_model_SUR_HC)
current_chrom = 2
plotEffectEstimate = function(current_chrom){
    hc_plot = ggplot(filter(effects, chrom == current_chrom), aes(marker, mean, group = trait)) +
        geom_point() + facet_grid(trait~type) +
        geom_hline(yintercept = 0) +
        geom_point(size = 0.3) +
        #geom_point(data = filter(true_effects, chrom == current_chrom), aes(y = true_effects), color = "red") +
        geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, size = 0.3)
    save_plot(paste0("data/figures/stan_SUR_HC_chrom", current_chrom, ".png"), hc_plot, base_height = 6, base_aspect_ratio = 1.8)
    return(hc_plot)
}
llply(seq_along(loci_per_chrom), plotEffectEstimate)
