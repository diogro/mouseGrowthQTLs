setwd("/home/diogro/projects/mouse-qtls")
source('read_mouse_data.R')
source('utils.R')

Rdatas_folder = "~/gdrive/LGSM_project_Rdatas/"
#Rdatas_folder = "./data/Rdatas/"

library(rstan)

area_data = inner_join(area_phen_std,
                       simulated_genomes[[2]],
                       by = "ID")

area_data = area_data %>% mutate_each(funs(scale), matches('area'))

true_effects = c(findA(area_data$area1, area_data$chrom4_A10, 0.05),
                 findA(area_data$area2, area_data$chrom4_A10, 0.05),
                 findA(area_data$area3, area_data$chrom4_A13, 0.05),
                 findA(area_data$area4, area_data$chrom4_A14, 0.05),
                 findA(area_data$area5, area_data$chrom4_A14, 0.05),
                 findA(area_data$area6, area_data$chrom4_A13, 0.05),
                 findA(area_data$area7, area_data$chrom4_A14, 0.05))
true_effects = data.frame(true_effects, trait = area_traits, type = "additive",
                          marker = c(rep(10, 2), 13, rep(14, 2), 13, 14))

area_data$area1 = makeSimData(area_data$area1, area_data$chrom4_A10, 0.05)
area_data$area2 = makeSimData(area_data$area2, area_data$chrom4_A10, 0.05)
area_data$area3 = makeSimData(area_data$area3, area_data$chrom4_A13, 0.05)
area_data$area4 = makeSimData(area_data$area4, area_data$chrom4_A14, 0.05)
area_data$area5 = makeSimData(area_data$area5, area_data$chrom4_A14, 0.05)
area_data$area6 = makeSimData(area_data$area6, area_data$chrom4_A13, 0.05)
area_data$area7 = makeSimData(area_data$area7, area_data$chrom4_A14, 0.05)

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
    param_list = list(K        =  K,
                      J        =  J,
                      N        =  N,
                      n_family =  n_family,
                      family   =  family,
                      ad       =  ad,
                      dm       =  dm,
                      y        =  y)
    return(param_list)
}
stan_model = stan(file = './mixedModelGmatrix.stan',
                  data = getStanInput(4), chain=4, iter = 1000)

pairs(stan_model)

getStanEffects = function(stan_model){
  HC_summary = summary(stan_model)
  s = loci_per_chrom[current_chrom] * num_area_traits * 2
  effects = data.frame(HC_summary$summary[1:s, c(1, 4, 7)])
  colnames(effects) <- c("mean", "lower", "upper")
  effects$type = rep(c("additive", "dominance"), each = s/2)
  effects$chrom = current_chrom
  effects$marker = rep(1:loci_per_chrom[current_chrom], 2*num_area_traits)
  effects$trait = rep(area_traits, each = loci_per_chrom[current_chrom])
  tbl_df(effects)
}
effects = getStanEffects(stan_model)
hc_plot = ggplot(effects, aes(marker, mean, group = trait)) +
  geom_point() + facet_grid(trait~type, scales = "free") +
  geom_hline(yintercept = 0) +
  geom_point(size = 0.3) +
  geom_point(data = true_effects, aes(y = true_effects), color = "red") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, size = 0.3)
save_plot("data/figures/sim_stan_sparceRegression.png", hc_plot, base_height = 6, base_aspect_ratio = 1.8)
