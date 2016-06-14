if(!require(install.load)) {install.packages("install.load"); library(install.load)}
install_load("plyr", "evolqg", "dplyr", "tidyr", "readr", "ggplot2", "cowplot")

getEffects = function(trait){
list(
  single_eff = read_csv(paste0("./data/", trait," traits/effectsSingleLocus.csv")),
  i5cM_eff   = read_csv(paste0("./data/", trait," traits/", trait,"_effectsInterval_5cM_mcmc.csv")),
  i10cM_eff  = read_csv(paste0("./data/", trait," traits/", trait,"_effectsInterval_10cM_mcmc.csv")),
  i15cM_eff  = read_csv(paste0("./data/", trait," traits/", trait,"_effectsInterval_15cM_mcmc.csv")),
  i20cM_eff  = read_csv(paste0("./data/", trait," traits/", trait,"_effectsInterval_20cM_mcmc.csv")),
  singleDIC  = read_csv(paste0("./data/", trait," traits/singleLocusDIC.csv")),
  intDIC     = read_csv(paste0("./data/", trait," traits/intervalMappingDIC.csv")) 
 # nullDIC    = readRDS(file = paste0("./data/", trait," traits/nullDIC.rds"))
)}
trait_sets = c("area", "growth", "necropsy")
effects_list = llply(trait_sets, getEffects)
names(effects_list) = trait_sets

effects_list[['necropsy']]$i10cM_eff %>% 
  select(count, chrom, marker, trait, ad_post.mean, dm_post.mean) %>%
  rename(additive = ad_post.mean, dominance = dm_post.mean) %>%
  gather(type, value, additive:dominance) %>% 
  mutate(trait = as.factor(trait)) %>%
  ggplot(aes(count, value, group = interaction(trait, type), color = trait)) +
  geom_hline(yintercept = 0) + 
  geom_line() + facet_grid(trait~type) 

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
