source('./read_mouse_data.R')
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
trait_sets = c("area", "growth", "necropsy", "weight")
effects_list = llply(trait_sets, getEffects)
names(effects_list) = trait_sets

chrom_limits = inner_join(select(effects_list[['necropsy']]$single, chrom, marker, count),
           effects_list[[1]]$singleDIC %>% group_by(chrom) %>% summarise(marker = max(locus)),
           by = c("chrom", "marker")) %>% mutate(count = as.numeric(count) + 0.5) %>% unique

current_chrom = 6
x = effects_list[[4]]$i10cM
current_trait = "weight"
plotMainEffects = function(x, current_chrom, title = NULL, alpha = TRUE){
  plot_data =
    x %>%
    select(count, chrom, marker, trait,
           ad_mean, ad_lower, ad_upper,
           dm_mean, dm_lower, dm_upper, DICDiff) %>%
    filter(chrom == current_chrom) %>%
    rename(additive = ad_mean, dominance = dm_mean) %>%
    gather(type, value, additive, dominance) %>%
    mutate(trait = as.factor(eval(parse(text = paste0(current_trait, "_traits")))[trait]),
           chrom = as.factor(chrom))

    plot_data %>%
        filter(chrom == current_chrom) %>%
        ggplot(aes(marker, value)) +
        geom_point(size = 2) +
        geom_errorbar(data = filter(plot_data, type == "additive"),
                      aes(ymax = ad_upper, ymin = ad_lower), width = 0) +
        geom_errorbar(data = filter(plot_data, type == "dominance"),
                      aes(ymax = dm_upper, ymin = dm_lower), width = 0) +
        geom_hline(yintercept = 0) +
        facet_grid(trait ~ type, scales = "free") +
        ggtitle(paste(title))
}
plotMainEffects(effects_list[[4]]$single, 3, alpha = FALSE)

plotAllEffects = function(current_trait, prefix = "", ...){
    effect_types = c("i5cM","i10cM","i15cM","i20cM", "single")
    llply(effect_types, function(x) {
              for(chrom in 1:19){
                  save_plot(paste0("data/figures/", prefix, current_trait, "_",x,"_chrom", chrom, ".png"),
                            plotMainEffects(effects_list[[current_trait]][[x]],
                                            chrom, paste(current_trait, x, "chrom", chrom), ...),
                            base_aspect_ratio = 2, base_height = 8)
              }
                      })
}

for(current_trait in trait_sets) plotAllEffects(current_trait, alpha = FALSE)
for(current_trait in trait_sets) plotAllEffects(current_trait, "alpha_", alpha = TRUE)


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
