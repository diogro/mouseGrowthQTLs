source("./read_mouse_data.R")
load("./data/Rdatas/all_loci_mcmc.Rdata")
all_effectsSingle_mcmc = read_csv("./data/area\ traits/effects_mcmc.csv")

install_load("evolqg")

dic_diff = (area_MCMC_null_model$DIC - unlist(lapply(all_loci_MCMC, function(x) x$DIC)))

sas_effects = read_csv("./data/area\ traits/Multivariate\ QTL\ analysi\ results/vectorsAreas_with\ sex_no\ i.csv")
sas_loci_test = read_csv("./data/area\ traits/Multivariate\ QTL\ analysi\ results/TESTSAreas_with\ sex_no\ i.csv")

p_ad = ggplot(data_frame(p = ifelse(all_effectsSingle_mcmc$p_ad > 0.5,
                                    all_effectsSingle_mcmc$p_ad,
                                    1 - all_effectsSingle_mcmc$p_ad)) %>%
                mutate(sas = sas_effects$Additive,
                       sas_upper = ifelse(sas_effects$aLOD > 2 | sas_effects$aLOD < -log10(0.99),
                                          sas_effects$Additive + 2 * sas_effects$AdditiveSE, NA),
                       sas_lower = ifelse(sas_effects$aLOD > 2 | sas_effects$aLOD < -log10(0.99),
                                          sas_effects$Additive - 2 * sas_effects$AdditiveSE, NA),
                       mcmc = all_effectsSingle_mcmc$ad_post.mean,
                       mcmc_lower = ifelse(p > 0.99, all_effectsSingle_mcmc$ad_l.95..CI, NA),
                       mcmc_upper = ifelse(p > 0.99, all_effectsSingle_mcmc$ad_u.95..CI, NA)),
              aes(sas, mcmc)) +
  geom_point(aes(alpha = p)) +
  geom_errorbarh(size = 0.3, color = "orange", aes(xmax = sas_upper, xmin = sas_lower, height = 0, alpha = p)) +
  geom_errorbar(size = 0.3, color = "blue", aes(ymax = mcmc_upper, ymin = mcmc_lower, height = 0, alpha = p)) +
  geom_smooth(method = "lm") +
  geom_abline() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
  scale_alpha_continuous(range = c(0, 1)) +
  ggtitle("Additive")

p_dm = ggplot(data_frame(p = ifelse(all_effectsSingle_mcmc$p_dm > 0.5,
                                    all_effectsSingle_mcmc$p_dm,
                                    1 - all_effectsSingle_mcmc$p_dm)) %>%
                mutate(sas = sas_effects$Dominance,
                       sas_upper = ifelse(sas_effects$dLOD > 2 | sas_effects$dLOD < -log10(0.99),
                                          sas_effects$Dominance + 2 * sas_effects$DominanceSE, NA),
                       sas_lower = ifelse(sas_effects$dLOD > 2 | sas_effects$dLOD < -log10(0.99),
                                          sas_effects$Dominance - 2 * sas_effects$DominanceSE, NA),
                       mcmc = all_effectsSingle_mcmc$dm_post.mean,
                       mcmc_lower = ifelse(p > 0.99, all_effectsSingle_mcmc$dm_l.95..CI, NA),
                       mcmc_upper = ifelse(p > 0.99, all_effectsSingle_mcmc$dm_u.95..CI, NA)),
              aes(sas, mcmc)) +
  geom_point(aes(alpha = p)) +
  geom_errorbarh(size = 0.3, color = "orange", aes(xmax = sas_upper, xmin = sas_lower, height = 0, alpha = p)) +
  geom_errorbar(size = 0.3, color = "blue", aes(ymax = mcmc_upper, ymin = mcmc_lower, height = 0, alpha = p)) +
  geom_smooth(method = "lm") +
  geom_abline() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
  scale_alpha_continuous(range = c(0, 1)) + ggtitle("Dominance")

p_aLOD = ggplot(data.frame(sas = sas_effects$aLOD, mcmc = -log10(1-all_effectsSingle_mcmc$p_ad)),
       aes(sas, mcmc)) + geom_point() + geom_abline() + ggtitle("Additive LPR-values")
p_dLOD = ggplot(data.frame(sas = sas_effects$dLOD, mcmc = -log10(1-all_effectsSingle_mcmc$p_dm)),
       aes(sas, mcmc)) + geom_point() + geom_abline() + ggtitle("Dominance LPR-values")
p_pLoci = ggplot(data.frame(sas = sas_loci_test$mLOD, mcmc = dic_diff, CHR = as.factor(sas_loci_test$CHR)),
                 aes(sas, mcmc, color = CHR)) + geom_text(aes(label = CHR)) +
                 geom_smooth(method = "lm", color = "blue") +
                 geom_vline(xintercept = 6) +
                 ggtitle("mLOD - DIC")


cf = coef(lm(dic_diff~sas_loci_test$mLOD))

locus_sig = gather(data.frame(locus = sas_loci_test$Count,
                              chrom = as.factor(sas_loci_test$CHR),
                              sas = sas_loci_test$mLOD,
                              dic = dic_diff), variable, value, sas:dic)
mLOD_plot = ggplot(locus_sig,
                   aes(locus, value, group = variable, color = chrom)) +
  geom_line() +
  geom_hline(data = data_frame(cutoff = c(cf[1] + 5*cf[2], 5), variable = c("dic", "sas")), aes(yintercept = cutoff)) +
  facet_wrap(~variable, ncol = 1, scales = "free_y") + ggtitle("Per marker test")

mcmc_sas_comparions_plot = plot_grid(p_ad, p_dm, p_pLoci, p_aLOD, p_dLOD)
save_plot("./data/area traits/mcmc_sas_comparison.png", mcmc_sas_comparions_plot, ncol = 3, nrow = 2, base_height = 5)
save_plot("./data/area traits/mcmc_sas_comparison_LOD.png", mLOD_plot, base_height = 7, base_aspect_ratio= 2)


