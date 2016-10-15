library(rstan)
rstan_options(auto_write = TRUE)
getStanInputFullGenome = function(current_data, trait_vector, J = sum(loci_per_chrom))
{
    K        = length(trait_vector)
    N        = dim(current_data)[1]
    n_family = length(unique(current_data$FAMILY))
    family   = as.integer(as.factor(current_data$FAMILY))
    ad       = as.matrix(select(current_data, matches(paste0('chrom.*_A'))))
    dm       = as.matrix(select(current_data, matches(paste0('chrom.*_D'))))
    teVec_ad = t(eigen(cov(ad))$vectors)
    teVec_dm = t(eigen(cov(dm))$vectors)
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
                      teVec_ad = teVec_ad,
                      teVec_dm = teVec_dm,
                      y        = y,
                      beta_ad  = beta_ad,
                      beta_dm  = beta_dm)
    return(param_list)
}
getStanInput = function(current_chrom, current_data, trait_vector,
                        J = loci_per_chrom[current_chrom])
{
    K        = length(trait_vector)
    N        = dim(current_data)[1]
    n_family = length(unique(current_data$FAMILY))
    family   = as.integer(as.factor(current_data$FAMILY))
    ad       = as.matrix(select(current_data, matches(paste0('chrom', current_chrom, '_A'))))
    dm       = as.matrix(select(current_data, matches(paste0('chrom', current_chrom, '_D'))))
    teVec_ad = t(eigen(cov(ad))$vectors)
    teVec_dm = t(eigen(cov(dm))$vectors)
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
                      teVec_ad = teVec_ad,
                      teVec_dm = teVec_dm,
                      y        = y,
                      beta_ad  = beta_ad,
                      beta_dm  = beta_dm)
    return(param_list)
}
getStanEffects = function(current_chrom, stan_model, trait_vector,
                          J = loci_per_chrom[current_chrom],
                          markers = 1:loci_per_chrom[current_chrom])
{
  K = length(trait_vector)
  HC_summary = summary(stan_model, pairs = c("w_ad", "w_dm"))$summary
  s = J * K * 2
  mask = grepl("w_", rownames(HC_summary))
  effects = data.frame(HC_summary[mask, c("mean", "2.5%", "97.5%")])
  colnames(effects) <- c("mean", "lower", "upper")
  effects$type = rep(c("additive", "dominance"), each = s/2)
  effects$chrom = current_chrom
  effects$marker = rep(markers, 2*K)
  effects$trait = rep(trait_vector, each = J)
  tbl_df(effects)
}
getStanShrinkage = function(current_chrom, stan_model, trait_vector,
                            J = loci_per_chrom[current_chrom],
                            markers = 1:loci_per_chrom[current_chrom])
{
  K = length(trait_vector)
  shrink = rstan::extract(stan_model, pars = c("shrink_ad", "shrink_dm"))
  raw_weights_ad = shrink[[1]]
  raw_weights_dm = shrink[[2]]
  weights_ad = reshape2::melt(raw_weights_ad) %>% ddply(.(Var2, Var3), numcolwise(mean))
  weights_dm = reshape2::melt(raw_weights_dm) %>% ddply(.(Var2, Var3), numcolwise(mean))
  weights = rbind(select(weights_ad, -iterations), select(weights_dm, -iterations))
  s = J * K * 2
  colnames(weights) <- c("trait", "marker", "mean")
  weights$type = rep(c("additive", "dominance"), each = s/2)
  weights$chrom = current_chrom
  weights$marker = rep(markers, 2*K)
  weights$trait = rep(trait_vector, each = J)
  tbl_df(weights)
}
plotEffectEstimate = function(current_chrom, effects, file = NULL, true_effects = NULL, scale = "fixed")
{
    hc_plot = ggplot(filter(effects, chrom == current_chrom), aes(marker, mean, group = trait)) +
        geom_point() + facet_grid(trait~type, scale = scale) +
        geom_hline(yintercept = 0) +
        geom_point(size = 0.3) + ggtitle(paste("Effects chrom", current_chrom)) +
        geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, size = 0.3)
    if(!is.null(true_effects))
        hc_plot = hc_plot + geom_point(data = filter(true_effects, chrom == current_chrom), aes(y = true_effects), color = "red")
    if(!is.null(file))
        save_plot(paste0("./data/figures/", file, "_effects_chrom", current_chrom, ".png"),
                  hc_plot, base_height = 6, base_aspect_ratio = 1.8)
    return(hc_plot)
}
plotShrinkage = function(current_chrom, weights, file = NULL)
{
    hc_plot = ggplot(filter(weights, chrom == current_chrom), aes(marker, mean, group = trait)) +
        geom_point() + facet_grid(trait~type, scale = "free") +
        geom_hline(yintercept = 0) +
        geom_hline(yintercept = 0.5, linetype = "dashed") +
        geom_point(size = 0.3) + ggtitle(paste("Shrinkage chrom", current_chrom))
    if(!is.null(file))
        save_plot(paste0("./data/figures/", file, "_shrinkage_chrom", current_chrom, ".png"),
                  hc_plot, base_height = 6, base_aspect_ratio = 1.8)
    return(hc_plot)
}
runStanModel = function(current_chrom, current_data, trait_vector, chain = 4, iter = 200, warmup = 100,
                        model_file = './SUR_horseShoePlus.stan', ...)
{
    stan_model_SUR_HC = stan(file = model_file,
                             data = getStanInput(current_chrom, current_data, trait_vector),
                             chain=chain, iter = iter, warmup = warmup, ...)
    effects = getStanEffects(current_chrom, stan_model_SUR_HC, trait_vector)
    weights = getStanShrinkage(current_chrom, stan_model_SUR_HC, trait_vector)
    return(list(effects, weights, rstan::extract(stan_model_SUR_HC)))
}
runStanModelFullGenome = function(current_data, trait_vector, iter = 200, warmup = 100, parallel = TRUE,
                                  model_file = "./SUR_horseShoePlus.stan", ...)
{
    all_chroms = alply(1:19, 1, runStanModel, current_data, trait_vector,
                       chain = 1, iter = iter, model_file = model_file, ..., .parallel = parallel, .inform = TRUE)
    effects = Reduce(bind_rows, llply(all_chroms, '[[', 1))
    weights = Reduce(bind_rows, llply(all_chroms, '[[', 2))
    a_fits = llply(all_chroms, '[[', 3)
    return(list(effects = effects, weights = weights, fits = a_fits))
}


