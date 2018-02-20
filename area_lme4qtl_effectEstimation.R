setwd("/home/diogro/projects/mouse-qtls")
source('read_mouse_data.R')
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#Rdatas_folder = "~/gdrive/LGSM_project_Rdatas/"
Rdatas_folder = "./Rdatas/"

area_data = inner_join(area_phen_std,
                        area_markers,
                        by = "ID") %>%
  gather(variable, value, area1:area7)

area_data = mutate(area_data, FAMILY = as.factor(FAMILY), variable = as.factor(variable))
null_formula = "value ~ variable + (0 + variable|FAMILY)"

makeMarkerList = function(pos) paste(paste('variable:chrom', pos[1],"_", c('A', 'D'), pos[2], sep = ''), collapse = ' + ')
markerMatrix = ldply(1:19, function(x) data.frame(chrom = x, marker = 1:loci_per_chrom[[x]]))
markerList = alply(markerMatrix, 1, makeMarkerList)
markerPositions = cbind(markerMatrix, read_csv("./data/markers/marker_positions.csv")[,3])
names(markerPositions)[3] = "cM"

significantMarkerMatrix = read_csv("./data/area_significant_markers.csv")

significantMarkerList = alply(significantMarkerMatrix, 1, makeMarkerList)
significant_marker_term = paste(significantMarkerList, collapse = " + ")

genotype_formula = paste(null_formula, significant_marker_term, sep = ' + ')

#significant_markers_model = lmer(as.formula(genotype_formula),
									   #data = area_data,
									   #REML = FALSE)
#save(significant_markers_model, file = "./Rdatas/area_significant_lmer_fit.Rdata")
load(file = "./Rdatas/area_significant_lmer_fit.Rdata")

coef_df = summary(significant_markers_model)$coef[-c(1:7), ]
coef_mean = matrix(coef_df[,1], ncol = 7, byrow = T)
coef_sd = matrix(coef_df[,2], ncol = 7, byrow = T)
ID = unlist(paste(rep(c("A", "D"), nrow(significantMarkerMatrix)),
                  rep(apply(significantMarkerMatrix, 1, paste, collapse = "_"), each = 2)) %>% {gsub(" ", "", .)})

colnames(coef_mean) = colnames(coef_sd) = area_traits
rownames(coef_mean) = rownames(coef_sd) = ID 

coef_upper = coef_mean + 2*coef_sd
coef_lower = coef_mean - 2*coef_sd

ID = unlist(rep(apply(significantMarkerMatrix, 1, paste, collapse = "_"), each = 2) %>% {gsub(" ", "", .)})

effects = data.frame(coef_mean) %>%
           mutate(id = ID) %>%
           gather(trait, mean, area1:area7)
effects$class = c("additive", "dominance")
effects$upper = (data.frame(coef_upper) %>% mutate(id = ID) %>% gather(trait, upper, area1:area7))$upper
effects$lower = (data.frame(coef_lower) %>% mutate(id = ID) %>% gather(trait, lower, area1:area7))$lower
write_csv(effects, "./data/area_significant_marker_effects.csv")

effects_plot = ggplot(effects, aes(id, mean)) + geom_point() + facet_wrap(~class) + geom_pointrange(aes(ymin = lower, ymax = upper)) + facet_grid(trait~class) + geom_hline(yintercept = 0)

getGenColumn <- function(pos, type){
    chrom = pos[1]
    marker = pos[2]
    gsub(" ", "", paste0("chrom", chrom, "_", type, marker))
}

area_data_wide = inner_join(area_phen_std,
                              area_markers,
                              by = "ID")
K        = length(area_traits)
N        = dim(area_data_wide)[1]
J        = nrow(significantMarkerMatrix)
n_family = length(unique(area_data_wide$FAMILY))
family   = as.integer(as.factor(area_data_wide$FAMILY))
ad       = as.matrix(dplyr::select(area_data_wide, apply(significantMarkerMatrix, 1, getGenColumn, "A")))
dm       = as.matrix(dplyr::select(area_data_wide, apply(significantMarkerMatrix, 1, getGenColumn, "D")))
y        = as.matrix(dplyr::select(area_data_wide, area_traits))
param_list = list(K        = K,
                  J        = J,
                  N        = N,
                  n_family = n_family,
                  family   = family,
                  ad       = ad,
                  dm       = dm,
                  y        = y)
                  
stan_model_SUR_HC = stan(file = "./SUR_horseShoePlus.stan",
                         data = param_list,
                         chain=6, iter = 700, warmup = 500,
                         control = list(adapt_delta = 0.99, max_tree_depth = 12))
save(stan_model_SUR_HC, file = "./Rdatas/area_significant_stan_HCp_fit.Rdata")
load(file = "./Rdatas/area_significant_stan_HCp_fit.Rdata")
plot(stan_model_SUR_HC, pars = "w_ad")

stan_model_SUR = stan(file = "./SUR.stan",
                         data = param_list,
                         chain=4, iter = 200, warmup = 100,
                         control = list(adapt_delta = 0.99, max_tree_depth = 12))

save(stan_model_SUR, file = "./Rdatas/area_significant_stan_fit.Rdata")
load(file = "./Rdatas/area_significant_stan_fit.Rdata")
plot(stan_model_SUR, pars = "w_ad")

getStanEffects = function(markerMatrix, stan_model, trait_vector,
                          J = nrow(markerMatrix),
                          markers = 1:loci_per_chrom[current_chrom])
{
  K = length(trait_vector)
  HC_summary = summary(stan_model, pairs = c("w_ad", "w_dm"))$summary
  s = J * K * 2
  mask = grepl("w_", rownames(HC_summary))
  effects = data.frame(HC_summary[mask, c("mean", "2.5%", "97.5%")])
  colnames(effects) <- c("mean", "lower", "upper")
  effects$class = rep(c("additive", "dominance"), each = s/2)
  effects$chrom = rep(markerMatrix$chrom, K)
  effects$marker = rep(markerMatrix$marker, 2*K)
  effects$trait = rep(trait_vector, each = J)
  effects$id = paste(effects$chrom, effects$marker, sep = "_")
  tbl_df(effects)
}

effectsStan = getStanEffects(significantMarkerMatrix, stan_model_SUR_HC, area_traits)
effectsStan_plot = ggplot(effectsStan, aes(id, mean)) + geom_point() + facet_wrap(~class) + geom_pointrange(aes(ymin = lower, ymax = upper)) + facet_grid(trait~class) + geom_hline(yintercept = 0)

names(effects)
names(effectsStan)

plot_grid(effects_plot, effectsStan_plot)
write_csv(effectsStan[,names(effects)], "./data/area_significant_marker_effects_SUR.csv")

