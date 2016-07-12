getStanInput = function(current_chrom, current_data, trait_vector)
{
    K        = length(trait_vector)
    J        = loci_per_chrom[current_chrom]
    N        = dim(current_data)[1]
    n_family = length(unique(current_data$FAMILY))
    family   = as.integer(as.factor(current_data$FAMILY))
    ad       = as.matrix(select(current_data, matches(paste0('chrom', current_chrom, '_A'))))
    dm       = as.matrix(select(current_data, matches(paste0('chrom', current_chrom, '_D'))))
    y        = as.matrix(current_data[trait_vector])
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
getStanEffects = function(current_chrom, stan_model, trait_vector)
{
  K = length(trait_vector)
  HC_summary = summary(stan_model, pairs = c("w_ad", "w_dm"))$summary
  s = loci_per_chrom[current_chrom] * K * 2
  mask = grepl("w_", rownames(HC_summary))
  effects = data.frame(HC_summary[mask, c("mean", "2.5%", "97.5%")])
  colnames(effects) <- c("mean", "lower", "upper")
  effects$type = rep(c("aditive", "dominance"), each = s/2)
  effects$chrom = current_chrom
  effects$marker = rep(1:loci_per_chrom[current_chrom], 2*K)
  effects$trait = rep(trait_vector, each = loci_per_chrom[current_chrom])
  tbl_df(effects)
}
getStanShrinkage = function(current_chrom, stan_model, trait_vector)
{
  K = length(trait_vector)
  shrink = (extract(stan_model, pars = c("lambda_ad", "lambda_dm", "tau")))
  raw_weights_ad = sweep(shrink[[1]], 1, shrink[[3]], "*")
  raw_weights_dm = sweep(shrink[[2]], 1, shrink[[3]], "*")
  weights_ad = reshape2::melt(raw_weights_ad) %>% ddply(.(Var2, Var3), numcolwise(mean))
  weights_dm = reshape2::melt(raw_weights_dm) %>% ddply(.(Var2, Var3), numcolwise(mean))
  weights = rbind(select(weights_ad, -iterations), select(weights_dm, -iterations))
  s = loci_per_chrom[current_chrom] * K * 2
  colnames(weights) <- c("trait", "marker", "mean")
  weights$type = rep(c("aditive", "dominance"), each = s/2)
  weights$chrom = current_chrom
  weights$marker = rep(1:loci_per_chrom[current_chrom], 2*K)
  weights$trait = rep(trait_vector, each = loci_per_chrom[current_chrom])
  tbl_df(weights)
}
plotEffectEstimate = function(current_chrom, effects, file = NULL)
{
    hc_plot = ggplot(filter(effects, chrom == current_chrom), aes(marker, mean, group = trait)) +
        geom_point() + facet_grid(trait~type, scales = "free") +
        geom_hline(yintercept = 0) +
        geom_point(size = 0.3) +
        geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, size = 0.3)
    if(!is.null(file))
        save_plot(file, hc_plot, base_height = 6, base_aspect_ratio = 1.8)
    return(hc_plot)
}
plotWeights = function(current_chrom, weights, file = NULL)
{
    hc_plot = ggplot(filter(weights, chrom == current_chrom), aes(marker, mean, group = trait)) +
        geom_point() + facet_grid(trait~type) +
        geom_hline(yintercept = 0) +
        geom_hline(yintercept = 0.5, linetype = "dashed") + 
        geom_point(size = 0.3)
    if(!is.null(file))
        save_plot(file, hc_plot, base_height = 6, base_aspect_ratio = 1.8)
    return(hc_plot)
}
runStanModel = function(current_chrom, current_data, trait_vector, chain = 1, iter = 500, ...)
{
    stan_model_SUR_HC = stan(file = './SUR_horseShoe.stan',
                             data = getStanInput(current_chrom, current_data, trait_vector),
                             chain=chain, iter = iter, ...)
    effects = getStanEffects(current_chrom, stan_model_SUR_HC, trait_vector)
    weights = getStanShrinkage(current_chrom, stan_model_SUR_HC, trait_vector)
    return(list(effects, weights, rstan::extract(stan_model_SUR_HC)))
}

runStanModelFullGenome = function(current_data, trait_vector, iter = 2000, parallel = TRUE)
{
    all_chroms = llply(1:19, runStanModel, current_data, trait_vector, chain = 1, iter = iter, .parallel = parallel)
    effects = Reduce(bind_rows, llply(all_chroms, '[[', 1))
    weights = Reduce(bind_rows, llply(all_chroms, '[[', 2))
    models = llply(all_chroms, '[[', 3)
    return(list(effects = effects, weights = weights, models = models))
}

setwd("/home/diogro/projects/mouse-qtls")
source('read_mouse_data.R')
source('utils.R')

#Rdatas_folder = "~/gdrive/LGSM_project_Rdatas/"
Rdatas_folder = "./data/Rdatas/"

library(rstan)
rstan_options(auto_write = TRUE)

area_data = inner_join(area_phen_std,
                       Reduce(inner_join, area_markers),
                       by = "ID")

area_data = area_data %>% mutate_each(funs(scale), matches('area'))
area_data[area_traits] = scale(area_data[area_traits])
area_data[area_traits] = scale(area_data[area_traits], scale = rep(sqrt(1/20), num_area_traits))

library(doMC)
registerDoMC(19)
area_mapping = runStanModelFullGenome(area_data, area_traits, 2000, TRUE)
