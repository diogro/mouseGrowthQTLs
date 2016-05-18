setwd("../")
source("./read_mouse_data.R")
source("./qtlMapping/MCMCglmmScript.R")

install_load("evolqg")

G_sas = matrix(NA, 7, 7)
G_sas[1,1] = 0.0028634425
G_sas[2,1] = 0.005922578
G_sas[2,2] = 0.0525224556
G_sas[3,1] = 0.0019050943
G_sas[3,2] = 0.018124448
G_sas[3,3] = 0.0266428736
G_sas[4,1] = 0.0102800522
G_sas[4,2] = 0.0558587579
G_sas[4,3] = 0.0438159861
G_sas[4,4] = 0.1195128941
G_sas[5,1] = 0.0082103824
G_sas[5,2] = 0.0470287859
G_sas[5,3] = 0.0341308769
G_sas[5,4] = 0.086046372
G_sas[5,5] = 0.0828664803
G_sas[6,1] = 0.0013903878
G_sas[6,2] = 0.006217101
G_sas[6,3] = 0.0028513912
G_sas[6,4] = 0.0088679786
G_sas[6,5] = 0.0070761072
G_sas[6,6] = 0.0018280909
G_sas[7,1] = 0.007000065
G_sas[7,2] = 0.0320632749
G_sas[7,3] = 0.0186835342
G_sas[7,4] = 0.0520281583
G_sas[7,5] = 0.0465501872
G_sas[7,6] = 0.0055093826
G_sas[7,7] = 0.0337549975

G_sas[upper.tri(G_sas)] <- t(G_sas)[upper.tri(G_sas)]
G_mcmc/G_sas
MatrixCompare(G_sas, G_mcmc)
cov2cor(G_mcmc)
cov2cor(G_sas)

(cov(area_phen_std %>% select(area1:area7)) - (G_mcmc + R_mcmc))

dic_diff = (area_MCMC_null_model$DIC - unlist(lapply(all_loci_MCMC, function(x) x$DIC)))/4.76


sas_effects = read_csv("./data/area\ traits/Multivariate\ QTL\ analysi\ results/vectorsAreas_with\ sex_no\ i.csv")
sas_loci_test = read_csv("./data/area\ traits/Multivariate\ QTL\ analysi\ results/TESTSAreas_with\ sex_no\ i.csv")
hist(abs(sas_effects$Additive - all_effects_mcmc$ad_post.mean ))
p_ad = ggplot(data.frame(sas = sas_effects$Additive, mcmc = all_effects_mcmc$ad_post.mean),
       aes(sas, mcmc)) + geom_point() + geom_smooth(method = "lm") + geom_abline() + ggtitle("Additive")
p_dm = ggplot(data.frame(sas = sas_effects$Dominance, mcmc = all_effects_mcmc$dm_post.mean),
       aes(sas, mcmc)) + geom_point() + geom_smooth(method = "lm") + geom_abline() + ggtitle("Dominance")
p_cov = ggplot(data.frame(sas = G_sas[lower.tri(G_sas, diag = TRUE)], mcmc = G_mcmc[lower.tri(G_mcmc, diag = TRUE)]),
       aes(sas, mcmc)) + geom_point() + geom_smooth(method = "lm") + geom_abline() + geom_abline(slope = 2, color = "red") + ggtitle("Covariances")
p_cor = ggplot(data.frame(sas = cov2cor(G_sas)[lower.tri(G_sas, diag = TRUE)], mcmc = cov2cor(G_mcmc)[lower.tri(G_mcmc, diag = TRUE)]),
       aes(sas, mcmc)) + geom_point() + geom_smooth(method = "lm") + geom_abline() + ggtitle("Correlations")
p_aLOD = ggplot(data.frame(sas = sas_effects$aLOD, mcmc = -log10(1-all_effects_mcmc$p_am)),
       aes(sas, mcmc)) + geom_point() + geom_abline() + ggtitle("P-values")
p_dLOD = ggplot(data.frame(sas = 10^(-sas_effects$dLOD), mcmc = 1-all_effects_mcmc$p_dm),
       aes(sas, mcmc)) + geom_point() + geom_abline() + ggtitle("P-values")
p_pLoci = ggplot(data.frame(sas = sas_loci_test$mLOD, mcmc = dic_diff),
       aes(sas, mcmc)) + geom_point() + geom_smooth(method = "lm") + ggtitle("mLOD - DIC") + geom_abline()

locus_sig = gather(data.frame(locus = all_effects_mcmc$count, sas = sas_loci_test$mLOD, dic = dic_diff), variable, value, sas:dic)
mLOD_plot = ggplot(locus_sig, aes(locus, value, group = variable, color = variable)) + geom_line()

mcmc_sas_comparions_plot = plot_grid(p_ad, p_dm, p_cov, p_cor, p_aLOD, p_pLoci)
save_plot("./data/area traits/mcmc_sas_comparison.png", mcmc_sas_comparions_plot, ncol = 3, nrow = 2, base_height = 5)
save_plot("./data/area traits/mcmc_sas_comparison_LOD.png", mLOD_plot, base_height = 7, base_aspect_ratio= 2)


