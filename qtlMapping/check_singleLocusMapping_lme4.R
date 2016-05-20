source("./read_mouse_data.R")
source("./qtlMapping/area_traits_lme4_model.R")
source("./qtlMapping/check_G.R")

install_load("evolqg")


sas_effects = read_csv("./data/area\ traits/Multivariate\ QTL\ analysi\ results/vectorsAreas_with\ sex_no\ i.csv")
sas_loci_test = read_csv("./data/area\ traits/Multivariate\ QTL\ analysi\ results/TESTSAreas_with\ sex_no\ i.csv")

lme4_loci_sig = ldply(all_loci, function(x) x$p.value)$V1

p_ad = ggplot(data.frame(sas = sas_effects$Additive, lme4 = all_effects$ad),
       aes(sas, lme4)) + geom_point() + geom_smooth(method = "lm") + geom_abline() + ggtitle("Additive")
p_dm = ggplot(data.frame(sas = sas_effects$Dominance, lme4 = all_effects$dm),
       aes(sas, lme4)) + geom_point() + geom_smooth(method = "lm") + geom_abline() + ggtitle("Dominance")
p_aSE = ggplot(data.frame(sas = sas_effects$AdditiveSE, lme4 = all_effects$ad_se),
       aes(sas, lme4)) + geom_point() + geom_abline() + ggtitle("Standard Errors")
p_dSE = ggplot(data.frame(sas = sas_effects$DominanceSE, lme4 = all_effects$dm_se),
       aes(sas, lme4)) + geom_point() + geom_abline() + ggtitle("Standard Errors")
p_pLoci = ggplot(data.frame(sas = sas_loci_test$mPROB, lme4 = lme4_loci_sig),
       aes(sas, lme4)) + geom_point() + geom_smooth(method = "lm") + ggtitle("mLOD - DIC") + geom_abline()

locus_sig = gather(data.frame(locus = all_effects$count, sas = sas_loci_test$mLOD, lme4 = lme4_loci_sig), variable, value, sas:lme4)
mLOD_plot = ggplot(locus_sig, aes(locus, value, group = variable, color = variable)) + geom_line()

lme4_sas_comparions_plot = plot_grid(p_ad, p_dm, p_aSE, p_dSE)
save_plot("./data/area traits/lme4_sas_comparison.png", lme4_sas_comparions_plot, ncol = 3, nrow = 2, base_height = 5)
save_plot("./data/area traits/lme4_sas_comparison_LOD.png", mLOD_plot, base_height = 7, base_aspect_ratio= 2)


