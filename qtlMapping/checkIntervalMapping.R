 
source("./read_mouse_data.R")
load("./data/Rdatas/intervalMapping_15cM_mcmc.Rdata")
all_effectsInterval_mcmc = read_csv("./data/area traits/effectsInterval_15cM_mcmc.csv")

install_load("evolqg")

dic_diff = unlist(lapply(intervalMapping_MCMC, function(x) x$flanking$DIC / x$focal$DIC))

sas_effects = read_csv("./data/area traits/Multivariate QTL analysi results/vectorsAreasCIM_15cM.csv")
sas_loci_test = read_csv("./data/area traits/Multivariate QTL analysi results/TESTSAreasCIM_15cM.csv")

p_ad = ggplot(data_frame(p = ifelse(all_effectsInterval_mcmc$p_ad > 0.5, 
                                    all_effectsInterval_mcmc$p_ad, 
                                    1 - all_effectsInterval_mcmc$p_ad)) %>% 
                mutate(sas = sas_effects$Additive, 
                       sas_upper = ifelse(sas_effects$aLOD > 2 | sas_effects$aLOD < -log10(0.99), 
                                          sas_effects$Additive + 2 * sas_effects$AdditiveSE, NA),
                       sas_lower = ifelse(sas_effects$aLOD > 2 | sas_effects$aLOD < -log10(0.99), 
                                          sas_effects$Additive - 2 * sas_effects$AdditiveSE, NA),
                       mcmc = all_effectsInterval_mcmc$ad_post.mean,
                       mcmc_lower = ifelse(p > 0.99, all_effectsInterval_mcmc$ad_l.95..CI, NA),
                       mcmc_upper = ifelse(p > 0.99, all_effectsInterval_mcmc$ad_u.95..CI, NA)),
              aes(sas, mcmc)) + 
  geom_point(aes(alpha = p)) + 
  geom_errorbarh(size = 0.3, color = "orange", aes(xmax = sas_upper, xmin = sas_lower, height = 0, alpha = p)) +
  geom_errorbar(size = 0.3, color = "blue", aes(ymax = mcmc_upper, ymin = mcmc_lower, height = 0, alpha = p)) +
  geom_smooth(method = "lm") + 
  geom_abline() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
  scale_alpha_continuous(range = c(0, 1)) +
  ggtitle("Additive")

all_effectsInterval_mcmc %>% arrange(desc(ad_post.mean))

p_dm = ggplot(data_frame(p = ifelse(all_effectsInterval_mcmc$p_dm > 0.5, 
                                    all_effectsInterval_mcmc$p_dm, 
                                    1 - all_effectsInterval_mcmc$p_dm)) %>% 
                mutate(sas = sas_effects$Dominance, 
                       sas_upper = ifelse(sas_effects$dLOD > 2 | sas_effects$dLOD < -log10(0.99), 
                                          sas_effects$Dominance + 2 * sas_effects$DominanceSE, NA),
                       sas_lower = ifelse(sas_effects$dLOD > 2 | sas_effects$dLOD < -log10(0.99), 
                                          sas_effects$Dominance - 2 * sas_effects$DominanceSE, NA),
                       mcmc = all_effectsInterval_mcmc$dm_post.mean,
                       mcmc_lower = ifelse(p > 0.99, all_effectsInterval_mcmc$dm_l.95..CI, NA),
                       mcmc_upper = ifelse(p > 0.99, all_effectsInterval_mcmc$dm_u.95..CI, NA)),
              aes(sas, mcmc)) + 
  geom_point(aes(alpha = p)) + 
  geom_errorbarh(size = 0.3, color = "orange", aes(xmax = sas_upper, xmin = sas_lower, height = 0, alpha = p)) +
  geom_errorbar(size = 0.3, color = "blue", aes(ymax = mcmc_upper, ymin = mcmc_lower, height = 0, alpha = p)) +
  geom_smooth(method = "lm") + 
  geom_abline() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
  scale_alpha_continuous(range = c(0, 1)) + ggtitle("Dominance")

p_aLOD = ggplot(data.frame(sas = sas_effects$aLOD, mcmc = -log10(1-all_effectsInterval_mcmc$p_ad)),
       aes(sas, mcmc)) + geom_point() + geom_abline() + ggtitle("Additive LPR-values")
p_dLOD = ggplot(data.frame(sas = sas_effects$dLOD, mcmc = -log10(1-all_effectsInterval_mcmc$p_dm)),
       aes(sas, mcmc)) + geom_point() + geom_abline() + ggtitle("Dominance LPR-values")
p_pLoci = ggplot(data.frame(sas = sas_loci_test$mLOD, mcmc = dic_diff),
       aes(sas, mcmc)) + geom_point() + geom_smooth(method = "lm") + ggtitle("mLOD - DIC") + geom_abline()

lm(dic_diff ~ sas_loci_test$mLOD)
cor(dic_diff, sas_loci_test$mLOD)

locus_sig = gather(data.frame(locus = all_effectsInterval_mcmc$count,
                              chrom = as.factor(all_effectsInterval_mcmc$chrom),
                              sas = sas_loci_test$mLOD,
                              dic = dic_diff), variable, value, sas:dic)
mLOD_plot = ggplot(locus_sig, 
                   aes(locus, value, group = variable, color = chrom)) + 
  geom_line() + 
  facet_wrap(~variable, ncol = 1, scales = "free_y")

mcmc_sas_comparions_plot = plot_grid(p_ad, p_dm, p_pLoci, p_aLOD, p_dLOD)
save_plot("./data/area traits/mcmc_sas_comp_interval.png", mcmc_sas_comparions_plot, ncol = 3, nrow = 2, base_height = 5, base_aspect_ratio = 2)
save_plot("./data/area traits/mcmc_sas_comp_LOD_interval.png", mLOD_plot, base_height = 7, base_aspect_ratio= 2)