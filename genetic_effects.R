if(!require(install.load)) {install.packages("install.load"); library(install.load)}
install_load("plyr", "evolqg", "dplyr", "tidyr", "readr", "ggplot2", "cowplot")

area_single_eff  = read_csv("./data/area traits/effectsSingleLocus.csv")
area_int_5cM_eff = read_csv("./data/area traits/effectsInterval_5cM_mcmc.csv")
area_int_10cM_eff = read_csv("./data/area traits/effectsInterval_10cM_mcmc.csv")
area_int_15cM_eff = read_csv("./data/area traits/effectsInterval_15cM_mcmc.csv")
area_int_20cM_eff = read_csv("./data/area traits/effectsInterval_20cM_mcmc.csv")
area_singleDIC = read_csv("./data/area traits/singleLocusDIC.csv") 
area_intDIC = read_csv("./data/area traits/intervalMappingDIC.csv") 
area_nullDIC = readRDS(file = "./data/area traits/nullDIC.rds")

growth_single_eff  = read_csv("./data/growth traits/effectsSingleLocus.csv")
growth_int_5cM_eff = read_csv("./data/growth traits/growth_effectsInterval_5cM_mcmc.csv")
growth_int_10cM_eff = read_csv("./data/growth traits/growth_effectsInterval_10cM_mcmc.csv")
growth_int_15cM_eff = read_csv("./data/growth traits/growth_effectsInterval_15cM_mcmc.csv")
growth_int_20cM_eff = read_csv("./data/growth traits/growth_effectsInterval_20cM_mcmc.csv")
growth_singleDIC = read_csv("./data/growth traits/singleLocusDIC.csv") 
growth_intDIC = read_csv("./data/growth traits/intervalMappingDIC.csv") 
growth_nullDIC = readRDS(file = "./data/growth traits/nullDIC.rds")
max(growth_nullDIC)

area_int_10cM_eff %>% 
  select(count, chrom, marker, trait, ad_post.mean, dm_post.mean) %>%
  rename(additive = ad_post.mean, dominance = dm_post.mean) %>%
  gather(type, value, additive:dominance) %>% 
  mutate(trait = as.factor(trait)) %>%
  ggplot(aes(count, value, group = interaction(trait, type), color = trait)) +
  geom_hline(yintercept = 0) + 
  geom_line() + facet_grid(trait~type) 

area_intDIC  %>%
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
  geom_line() + geom_hline(yintercept = 20) + geom_hline(yintercept = 0) + facet_wrap(~interval, ncol = 1)
