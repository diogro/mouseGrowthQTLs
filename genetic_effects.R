if(!require(install.load)) {install.packages("install.load"); library(install.load)}
install_load("plyr", "evolqg", "dplyr", "tidyr", "readr", "ggplot2", "cowplot")

getEffects = function(trait){
  single_eff = read_csv(paste0("./data/", trait," traits/effectsSingleLocus.csv"))
  i5cM_eff   = read_csv(paste0("./data/", trait," traits/", trait,"_effectsInterval_5cM_mcmc.csv"))
  i10cM_eff  = read_csv(paste0("./data/", trait," traits/", trait,"_effectsInterval_10cM_mcmc.csv"))
  i15cM_eff  = read_csv(paste0("./data/", trait," traits/", trait,"_effectsInterval_15cM_mcmc.csv"))
  i20cM_eff  = read_csv(paste0("./data/", trait," traits/", trait,"_effectsInterval_20cM_mcmc.csv"))
  singleDIC  = read_csv(paste0("./data/", trait," traits/singleLocusDIC.csv"))
  intDIC     = read_csv(paste0("./data/", trait," traits/intervalMappingDIC.csv")) 
 # nullDIC    = readRDS(file = paste0("./data/", trait," traits/nullDIC.rds"))
  single_eff = left_join(single_eff, 
                         select(singleDIC, chrom, locus, singleLocus_DICDiff) %>%
                           rename(marker = locus),
                       by = c("chrom", "marker")) %>%
    rename(DICDiff  = singleLocus_DICDiff, 
           ad_mean  = ad_post.mean, ad_lower = ad_l.95..CI, ad_upper = ad_u.95..CI,
           dm_mean  = dm_post.mean, dm_lower = dm_l.95..CI, dm_upper = dm_u.95..CI)
  i5cM_eff = left_join(i5cM_eff, select(intDIC, chrom, marker, DICDiff_5cM),
           by = c("chrom", "marker")) %>%
   rename(DICDiff  = DICDiff_5cM, 
          ad_mean  = ad_post.mean, ad_lower = ad_l.95..CI, ad_upper = ad_u.95..CI,
          dm_mean  = dm_post.mean, dm_lower = dm_l.95..CI, dm_upper = dm_u.95..CI)
  i10cM_eff = left_join(i10cM_eff, select(intDIC, chrom, marker, DICDiff_10cM),
                       by = c("chrom", "marker")) %>%
    rename(DICDiff  = DICDiff_10cM, 
           ad_mean  = ad_post.mean, ad_lower = ad_l.95..CI, ad_upper = ad_u.95..CI,
           dm_mean  = dm_post.mean, dm_lower = dm_l.95..CI, dm_upper = dm_u.95..CI)
  i15cM_eff = left_join(i15cM_eff, select(intDIC, chrom, marker, DICDiff_15cM),
                       by = c("chrom", "marker")) %>%
    rename(DICDiff  = DICDiff_15cM, 
           ad_mean  = ad_post.mean, ad_lower = ad_l.95..CI, ad_upper = ad_u.95..CI,
           dm_mean  = dm_post.mean, dm_lower = dm_l.95..CI, dm_upper = dm_u.95..CI)
  i20cM_eff = left_join(i20cM_eff, select(intDIC, chrom, marker, DICDiff_20cM),
                       by = c("chrom", "marker")) %>%
    rename(DICDiff  = DICDiff_20cM, 
           ad_mean  = ad_post.mean, ad_lower = ad_l.95..CI, ad_upper = ad_u.95..CI,
           dm_mean  = dm_post.mean, dm_lower = dm_l.95..CI, dm_upper = dm_u.95..CI)
  list(single = single_eff, i5cM = i5cM_eff, i10cM = i10cM_eff, i15cM = i15cM_eff, 
       i20cM = i20cM_eff, singleDIC = singleDIC ,intDIC = intDIC
       #,nullDIC = nullDIC
       )
 }
trait_sets = c("area", "growth", "necropsy")
effects_list = llply(trait_sets, getEffects)
names(effects_list) = trait_sets

chrom_limits = inner_join(select(effects_list[['necropsy']]$single, chrom, marker, count),
           singleDIC %>% group_by(chrom) %>% summarise(marker = max(locus)),
           by = c("chrom", "marker")) %>% mutate(count = as.numeric(count) + 0.5) %>% unique

plotMainEffects = function(x, title = NULL, alpha = TRUE){
  plot_data =
    x %>% 
    select(count, chrom, marker, trait, 
           ad_mean, ad_lower, ad_upper,
           dm_mean, dm_lower, dm_upper, DICDiff) %>%
    filter(DICDiff > 0) %>%
    rename(additive = ad_mean, dominance = dm_mean) %>%
    gather(type, value, additive, dominance) %>% 
    mutate(trait = as.factor(eval(parse(text = paste0(current.trait, "_traits")))[trait]),
           chrom = as.factor(chrom))
  
    if(alpha){
      p = ggplot(plot_data, aes(count, value, group = interaction(trait, type), color = chrom, alpha = DICDiff))
    } else{ p = ggplot(plot_data, aes(count, value, group = interaction(trait, type), color = chrom)) }
    p +
    geom_hline(yintercept = 0) + 
    geom_point(size = 0.3) + 
    geom_errorbar(data = filter(plot_data, type == "additive"),
                 aes(ymin = ad_lower, ymax = ad_upper), width = 0, size = 0.3) +
    geom_errorbar(data = filter(plot_data, type == "dominance"),
                  aes(ymin = dm_lower, ymax = dm_upper), width = 0, size = 0.3) +
    geom_vline(data = chrom_limits, aes(xintercept = count), size = 0.2) +
    #geom_text(data = chrom_limits, aes(count - 2, y = 0.4, label = chrom)) +
    xlim(xmin = 1, xmax = 353) +
    facet_grid(trait~type, scales = "free") + ggtitle(title)
}
plotMainEffects(effects_list[[3]]$single, alpha = FALSE)

plotAllEffects = function(current.trait, prefix = "", ...){
  effect_types = c("i5cM","i10cM","i15cM","i20cM", "single")
  llply(effect_types, function(x) save_plot(paste0("data/figures/", prefix, current.trait, "_",x, ".png"), 
                                            plotMainEffects(effects_list[[current.trait]][[x]], 
                                                            paste(current.trait, x), ...), 
                                            base_aspect_ratio = 2, base_height = 8))
}

for(current.trait in trait_sets) plotAllEffects(current.trait, alpha = FALSE)
for(current.trait in trait_sets) plotAllEffects(current.trait, "alpha_", alpha = TRUE)


effects_list[['necropsy']]$intDIC  %>%
  select(chrom, marker, contains("DICDiff")) %>%
  mutate(count = seq(353)) %>% 
  rename( i5cM = DICDiff_5cM, 
         i10cM = DICDiff_10cM,
         i15cM = DICDiff_15cM,
         i20cM = DICDiff_20cM) %>% 
  gather(interval, value, i20cM:i5cM) %>%
  mutate(interval = factor(interval, levels = c("i5cM", "i10cM", "i15cM", "i20cM")),
         chrom = factor(chrom)) %>%
  ggplot(aes(count, value, group = interval, color = chrom)) +
  geom_line() + geom_hline(yintercept = 20) + geom_hline(yintercept = 0) + 
  facet_wrap(~interval, ncol = 1, scales = "free")
