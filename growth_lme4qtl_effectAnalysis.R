#setwd("/home/diogro/projects/mouse-qtls")
source('read_mouse_data.R')

ncores = 2
registerDoMC(ncores)
options(mc.cores = ncores)
setMKLthreads(ncores)

#Rdatas_folder = "~/gdrive/LGSM_project_Rdatas/"
Rdatas_folder = "./Rdatas/"

vectorCor = function(x, y) Normalize(x) %*% Normalize(y)
markerCov = function(marker1, marker2){
  marker1_col = growth_markers[,makeMarkerList(marker1)]
  marker2_col = growth_markers[,makeMarkerList(marker2)]
  cov(marker1_col, marker2_col)[1]
}
makeMarkerList = function(pos) paste('chrom', pos[1],"_", 'A', pos[2], sep = '')
lt = function(x, diag = TRUE) x[lower.tri(x, diag = diag)]
CalcInt = function(x) sd(eigen(cov2cor(x))$values)/(nrow(x) - 1)
markerMatrix = ldply(1:19, function(x) data.frame(chrom = x, marker = 1:loci_per_chrom[[x]]))
significantMarkerMatrix = read_csv("./data/growth_significant_markers.csv")

growth_data = inner_join(growth_phen_std,
                        growth_markers,
                        by = "ID") %>%
  gather(variable, value, growth12:growth78)
growth_data = mutate(growth_data, FAMILY = as.factor(FAMILY), variable = as.factor(variable))

e = read_csv("./data/growth_significant_marker_effects_SUR.csv") %>% arrange(trait, id, class)
load(file = "./Rdatas/significant_stan_fit.Rdata")
full_HCp = readRDS("./Rdatas/growth_scaled_allmarkers_HCPlus")

a_effect_matrix = e %>%
    select(id, class, trait, mean) %>%
    spread(trait, mean) %>%
    filter(class == "additive") %>% select(-class)
d_effect_matrix = e %>%
  select(id, class, trait, mean) %>%
  spread(trait, mean) %>%
  filter(class == "dominance") %>% select(-class)

w_ad_HC = rstan::extract(full_HCp, pars = "w_ad")[[1]]
w_dm_HC = rstan::extract(full_HCp, pars = "w_dm")[[1]]
effect_matrix_ad_HC = aaply(w_ad_HC, c(2, 3), mean)
effect_matrix_dm_HC = aaply(w_dm_HC, c(2, 3), mean)

a_effect_matrix_HC = data.frame(id = paste(markerMatrix$chrom, markerMatrix$marker, sep = "_"), 
                                as.data.frame(t(effect_matrix_ad_HC)))
colnames(a_effect_matrix_HC)[2:8] = growth_traits
d_effect_matrix_HC = data.frame(id = paste(markerMatrix$chrom, markerMatrix$marker, sep = "_"), 
                                as.data.frame(t(effect_matrix_dm_HC)))
colnames(d_effect_matrix_HC)[2:8] = growth_traits
 
PCbiplot <- function(PC, ids, x="PC1", y="PC2") {
  # PC being a prcomp object
  data <- data.frame(obsnames=ids, PC$x)
  plot <- ggplot(data, aes_string(x=x, y=y)) + geom_text(alpha=.4, size=3, aes(label=obsnames))
  plot <- plot + geom_hline(yintercept = 0, size=.2) + geom_vline(xintercept = 0, size=.2)
  datapc <- data.frame(varnames=rownames(PC$rotation), PC$rotation)
  mult <- min(
    (max(data[,y]) - min(data[,y])/(max(datapc[,y])-min(datapc[,y]))),
    (max(data[,x]) - min(data[,x])/(max(datapc[,x])-min(datapc[,x])))
  )
  datapc <- transform(datapc,
                      v1 = .7 * mult * (get(x)),
                      v2 = .7 * mult * (get(y))
  )
  plot <- plot + coord_equal() + geom_text(data=datapc, aes(x=v1, y=v2, label=varnames), size = 5, vjust=1, color="red")
  plot <- plot + geom_segment(data=datapc, aes(x=0, y=0, xend=v1, yend=v2), arrow=arrow(length=unit(0.2,"cm")), alpha=0.75, color="red")
  plot + theme_cowplot()
}
a_biplot = PCbiplot(prcomp(a_effect_matrix[,growth_traits]), ids = a_effect_matrix$id)
d_biplot = PCbiplot(prcomp(d_effect_matrix[,growth_traits]), ids = a_effect_matrix$id)
ad_biplot = plot_grid(a_biplot, d_biplot, labels = c("C", "D"))

a_biplot_HC = PCbiplot(prcomp(a_effect_matrix_HC[,growth_traits]), ids = a_effect_matrix_HC$id)
d_biplot_HC = PCbiplot(prcomp(d_effect_matrix_HC[,growth_traits]), ids = d_effect_matrix_HC$id)
ad_biplot_HC = plot_grid(a_biplot_HC, d_biplot_HC, labels = c("C", "D"))

LG = c(3.785,4.435,8.43,7.395,2.995,1.85,2.085)
SM = c(3.31 ,2.98,3.82,2.175,0.765,1.165,0.51)
F3 = sapply(growth_phen[,growth_traits], mean)

d_z = LG - SM
load("./Rdatas/growth_CovMatrices.Rdata")
growth_sds = apply(growth_phen[,growth_traits], 2, sd)
G = G_stan #* outer(growth_sds, growth_sds)
png("~/G_LGSM", width = 600, height = 600)
corrplot.mixed(cov2cor(G), upper = "ellipse")
dev.off()
plot(eigen(G)$values)
G_ext4 = ExtendMatrix(G, ret.dim = 4)[[1]]
G_ext5 = ExtendMatrix(G, ret.dim = 5)[[1]]
solve(G, d_z)
(beta = solve(G_ext4, d_z))
solve(G_ext5, d_z)

mean_a = colMeans(a_effect_matrix[,growth_traits])
mean_d = colMeans(d_effect_matrix[,growth_traits])
vectorCor(mean_a, d_z)
vectorCor(mean_d, d_z)
vectorCor(mean_a, beta)
vectorCor(mean_d, beta)

calcVa = function(i, a_effects, d_effects, markerMatrix){
  trait_sd = sapply(growth_phen[,growth_traits], sd)
  a_effects = a_effects * trait_sd
  d_effects = d_effects * trait_sd
  current_chrom = markerMatrix[i,1]
  focal_marker = markerMatrix[i,]
  focal_marker_col = growth_markers[,makeMarkerList(focal_marker)]
  n = dim(focal_marker_col)[1]
  genotype_freq = table(focal_marker_col)/n
  # Variance due to focal marker
  q = genotype_freq[1] + 1/2 * genotype_freq[2]
  p = genotype_freq[3] + 1/2 * genotype_freq[2]
  # additive contribution to Va
  V_a = 2*p*q * outer(a_effects[,i], a_effects[,i]) 
  # Dominance contribution to Va
  V_a = V_a + 2*p*q * (q - p)^2 * outer(d_effects[,i], d_effects[,i])
  # Additive by dominance contribution to Va
  V_a = V_a + 2*p*q * (q - p)  * (outer(a_effects[,i], d_effects[,i]) + 
                                  outer(d_effects[,i], a_effects[,i]))
  # Variance due to LD with focal marker
  for(j in 1:nrow(markerMatrix)){
    if (i != j) V_a = V_a + markerCov(focal_marker, markerMatrix[j,]) * 
                                      outer(a_effects[,i], a_effects[,j])
  }
  V_a
}
calcVd = function(i, d_effects, markerMatrix){
  trait_sd = sapply(growth_phen[,growth_traits], sd)
  d_effects = d_effects * trait_sd
  current_chrom = markerMatrix[i,1]
  focal_marker = markerMatrix[i,]
  focal_marker_col = growth_markers[,makeMarkerList(focal_marker)]
  n = dim(focal_marker_col)[1]
  genotype_freq = table(focal_marker_col)/n
  # Variance due to focal marker
  q = genotype_freq[1] + 1/2 * genotype_freq[2]
  p = genotype_freq[3] + 1/2 * genotype_freq[2]
  # additive contribution to Va
  V_d = (2*p*q)^2 * outer(d_effects[,i], d_effects[,i]) 
  # Variance due to LD with focal marker
  for(j in 1:nrow(markerMatrix)){
    d2ij = markerCov(focal_marker, markerMatrix[j,])
    if (i != j) V_d = V_d + d2ij^2 * outer(d_effects[,i], d_effects[,j])
  }
  V_d
}

w_ad = rstan::extract(stan_model_SUR, pars = "w_ad")[[1]]
effect_matrix_additive = aaply(w_ad, c(2, 3), mean)
w_dm = rstan::extract(stan_model_SUR, pars = "w_dm")[[1]]
effect_matrix_dominance = aaply(w_dm, c(2, 3), mean)
save(w_ad, w_dm, effect_matrix_additive, effect_matrix_dominance, file = "Rdatas/growth_add_dom_effectsMatrix.Rdata")
load(file = "Rdatas/growth_add_dom_effectsMatrix.Rdata")

#Va = laply(seq_along(1:400), function(i) colSums(laply(1:nrow(significantMarkerMatrix), calcVa, w_ad[i,,], w_dm[i,,], significantMarkerMatrix)), .parallel = TRUE)
#Vd = laply(seq_along(1:400), function(i) colSums(laply(1:nrow(significantMarkerMatrix), calcVd, w_dm[i,,], significantMarkerMatrix)), .parallel = TRUE)
#save(Va, Vd, file = paste0(Rdatas_folder, "VaVd_QTL.Rdata"))
load(paste0(Rdatas_folder, "VaVd_QTL.Rdata"))

Va_mean = aaply(Va, c(2, 3), mean)
Va_upper = aaply(Va, c(2, 3), quantile, 0.975)
Va_lower = aaply(Va, c(2, 3), quantile, 0.025)

Vd_mean  = aaply(Vd, c(2, 3), mean)
Vd_lower = aaply(Vd, c(2, 3), quantile, 0.025)
Vd_upper = aaply(Vd, c(2, 3), quantile, 0.975)

G = aaply(Gs_stan, c(2, 3), mean)
dimnames(G) = list(1:7, 1:7)
G_lower = aaply(Gs_stan, c(2, 3), quantile, 0.025)
G_upper = aaply(Gs_stan, c(2, 3), quantile, 0.975)

G_dam_lower = aaply(Gs_dam, c(2, 3), quantile, 0.025)
G_dam_upper = aaply(Gs_dam, c(2, 3), quantile, 0.975)

Vg = 1/2 * Va + 1/4 * Vd

Vg_mean  = aaply(Vg, c(2, 3), mean)
Vg_lower = aaply(Vg, c(2, 3), quantile, 0.025)
Vg_upper = aaply(Vg, c(2, 3), quantile, 0.975)

matrices <- list(FullSib = G_stan,
                 "Va QTL" = Va_mean,
                 "Vd QTL" = Vd_mean,
                 "1/2 Va + 1/4 Vd" = Vg_mean,
                 "FullSib cross-foster" = G_cf,
                 "FullSib non-cross-foster" = G_ncf,
                 #"G_dam" = G_dam,
                 #"G_nurse" = G_nurse,
                 "beta" = beta,
                 "delta Z" = d_z)
PrintMatrix(matrices)

png("./data/growth_family_Vg_FullSibG_correlation.png", width = 1500, height = 1500)
old.par = par()
par(mfrow = c(2, 2), cex=2, oma = c(0, 0, 0, 0))
corrplot.mixed(cov2cor(G),       upper = "ellipse", mar = c(0, 0, 1, 0))
corrplot.mixed(cov2cor(Vg_mean), upper = "ellipse", mar = c(0, 0, 1, 0))
corrplot.mixed(cov2cor(Va_mean), upper = "ellipse", mar = c(0, 0, 1, 0))
corrplot.mixed(cov2cor(Vd_mean), upper = "ellipse", mar = c(0, 0, 1, 0))
dev.off()
par(old.par)

write.csv(MatrixCompare(Vd_mean, G), file = "./data/TalkStuff/Vd_FamilyG_comparison.csv")
write.csv(MatrixCompare(Va_mean, G), file = "./data/TalkStuff/Va_FamilyG_comparison.csv")
write.csv(MatrixCompare(Vg_mean, G), file = "./data/TalkStuff/Vg_FamilyG_comparison.csv")
write.csv(MatrixCompare(Vg_mean, Vd_mean), file = "./data/TalkStuff/Va_Vd_comparison.csv")

write.csv(MatrixCompare(Vd_mean, G_dam), file = "./data/TalkStuff/Vd_DamG_comparison.csv")
write.csv(MatrixCompare(Va_mean, G_dam), file = "./data/TalkStuff/Va_DamG_comparison.csv")
write.csv(MatrixCompare(Vg_mean, G_dam), file = "./data/TalkStuff/Vg_DamG_comparison.csv")
write.csv(MatrixCompare(G, G_dam), file = "./data/TalkStuff/GFamily_GDam_comparison.csv")

write.csv(MatrixCompare(Vd_mean, G_nurse), file = "./data/TalkStuff/Vd_nurseG_comparison.csv")
write.csv(MatrixCompare(Va_mean, G_nurse), file = "./data/TalkStuff/Va_nurseG_comparison.csv")
write.csv(MatrixCompare(Vg_mean, G_nurse), file = "./data/TalkStuff/Vg_nurseG_comparison.csv")
write.csv(MatrixCompare(G, G_nurse), file = "./data/TalkStuff/GFamily_Gnurse_comparison.csv")
write.csv(MatrixCompare(G_dam, G_nurse), file = "./data/TalkStuff/Gdam_Gnurse_comparison.csv")

data.frame(Vg = diag(Vg_mean), G = diag(G)) %>% gather %>%
  ggplot(aes(c(1:7, 1:7), value, group = key, color = key)) + geom_point() + geom_line()
png("./data/growth_family_qtl_cov.png", width = 1500, height = 800)
par(mfrow = c(1, 1), cex=2)
plot(lt(G)~lt(Vg_mean), pch = 19, 
     ylab = "Family G-matrix covariances", xlab = "Genetic covariances predicted from QTLs (1/2 * Va + 1/4 * Vd)", 
     main = "", xlim = c(-0.03, 0.17), ylim = c(-0.22, 0.6))
segments(x0 = lt(Vg_lower), y0 = lt(G), x1 = lt(Vg_upper), y1 = lt(G))
segments(x0 = lt(Vg_mean), y0 = lt(G_lower), x1 = lt(Vg_mean), y1 = lt(G_upper))
points(diag(G)~diag(Vg_mean), col = "tomato3", pch = 19)
abline(lm(lt(G)~lt(Vg_mean)))
abline(0, 1, col = "blue")
text(0.15, 0.12, "Identity", col = "blue")
text(0.11, 0.35, "Variances", col = "tomato3")
text(0.05, -0.05, "Co-variances")
dev.off()

png("./data/growth_dam_qtl_cov.png", width = 1500, height = 800)
par(mfrow = c(1, 1), cex=2)
plot(lt(G_dam)~lt(Vg_mean), pch = 19, 
     ylab = "Dam G-matrix covariances", xlab = "Genetic covariances predicted from QTLs (1/2 * Va + 1/4 * Vd)", 
     main = "Growth traits", xlim = c(-0.03, 0.17), ylim = c(-0.22, 0.6))
segments(x0 = lt(Vg_lower), y0 = lt(G_dam), x1 = lt(Vg_upper), y1 = lt(G_dam))
segments(x0 = lt(Vg_mean), y0 = lt(G_dam_lower), x1 = lt(Vg_mean), y1 = lt(G_dam_upper))
points(diag(G_dam)~diag(Vg_mean), col = "tomato3", pch = 19)
abline(lm(lt(G_dam)~lt(Vg_mean)))
abline(0, 1, col = "blue")
text(0.15, 0.12, "Identity", col = "blue")
text(0.11, 0.35, "Variances", col = "tomato3")
text(0.05, -0.05, "Co-variances")
dev.off()


summary(lm(lt(G)~lt(Vg_mean)))
markerMatrix = ldply(1:19, function(x) data.frame(chrom = x, marker = 1:loci_per_chrom[[x]]))
Va_HC = colSums(laply(1:353, calcVa, effect_matrix_ad_HC, effect_matrix_dm_HC, markerMatrix))
old.par = par()
par(mfrow = c(1, 2), cex=2)
corrplot.mixed(cov2cor(G), upper = "ellipse", main = "\n\n\nFull-Sib G-matrix Correlation")
corrplot.mixed(cov2cor(Va_HC), upper = "ellipse", main = "\n\n\nVa_HC")
par(old.par)

library(viridis)

vectorCor(d_z, beta)
random_vec = matrix(rnorm(7*1000), 1000, 7)
quantile(abs(apply(random_vec, 1, vectorCor, rep(1, 7))), 0.95)
crss = data.frame(beta = apply(a_effect_matrix[,growth_traits], 1, vectorCor, beta),
                    dz = apply(a_effect_matrix[,growth_traits], 1, vectorCor,  d_z)) %>% gather(class, value, beta:dz)

a_corrs = data.frame(marker = 1:32, 
                     betaCorr = abs(apply(a_effect_matrix[,growth_traits], 1, vectorCor, beta)),
                     dzCorr = abs(apply(a_effect_matrix[,growth_traits], 1, vectorCor, d_z)),
                     norm = apply(a_effect_matrix[,growth_traits], 1, Norm))
d_corrs = data.frame(marker = 1:32,
                     betaCorr = abs(apply(d_effect_matrix[,growth_traits], 1, vectorCor, beta)),
                     dzCorr = abs(apply(d_effect_matrix[,growth_traits], 1, vectorCor, d_z)),
                     norm = apply(d_effect_matrix[,growth_traits], 1, Norm))
write.csv(a_corrs, "./data/growth_additive_correlations_beta_dZ.csv")
write.csv(d_corrs, "./data/growth_dominance_correlations_beta_dZ.csv")
ggplot(crss, aes(class, value, fill = class)) + geom_violin()
additive_beta = ggplot(a_corrs, aes(norm, betaCorr)) + geom_point() + geom_smooth(method = "lm", color = "black") + 
  labs(x = "Additive effect vector norm", y = expression(paste("Alligment with ", beta)))
dominance_beta = ggplot(d_corrs, aes(norm, betaCorr)) + geom_point() + geom_smooth(method = "lm", color = "black") + 
  labs(x = "Dominance effect vector norm", y = expression(paste("Alligment with ", beta)))
additive_dz = ggplot(a_corrs, aes(norm, dzCorr)) + geom_point() + geom_smooth(method = "lm", color = "black") + 
  labs(x = "Additive effect vector norm", y = expression("Alligment with divergence"))
dominance_dz = ggplot(d_corrs, aes(norm, dzCorr)) + geom_point() + geom_smooth(method = "lm", color = "black") + 
  labs(x = "Dominance effect vector norm", y = expression("Alligment with divergence"))
regressions = plot_grid(additive_beta, dominance_beta, additive_dz, dominance_dz, labels = LETTERS[1:4])
save_plot("data/growth_effect_aligment_regressions.png", regressions, base_height = 4, base_aspect_ratio = 2, ncol = 2, nrow = 2)

lm(beta~norm, data = a_corrs) %>% summary
lm(norm~d_z, data = a_corrs) %>% summary

growth_m = as.numeric(ddply(growth_phen, .(SEX), numcolwise(mean))[2,growth_traits])
growth_f = as.numeric(ddply(growth_phen, .(SEX), numcolwise(mean))[1,growth_traits])

posterior_predict = function(beta_ad){
    SM_e = rowSums(-1 * beta_ad) * sapply(growth_phen[,growth_traits], sd) + 
        sapply(growth_phen[,growth_traits], mean)
    LG_e = rowSums(beta_ad) * sapply(growth_phen[,growth_traits], sd) + 
        sapply(growth_phen[,growth_traits], mean)
    reshape2::melt(data.frame(trait = as.factor(growth_traits), 
                                              SM_Predicted = SM_e,
                                              LG_Predicted = LG_e)) %>% separate(variable, c("Line", "Type"))
}

w_ad = rstan::extract(stan_model_SUR, pars = "w_ad")[[1]]
effect_matrix = aaply(w_ad, c(2, 3), mean)
SM_e = rowSums(-1 * as.matrix(effect_matrix)) * sapply(growth_phen[,growth_traits], sd) + 
    sapply(growth_phen[,growth_traits], mean)
LG_e = rowSums(as.matrix(effect_matrix)) * sapply(growth_phen[,growth_traits], sd) + 
    sapply(growth_phen[,growth_traits], mean)
post = adply(w_ad, 1, posterior_predict)
growth_prediction = reshape2::melt(data.frame(trait = as.factor(growth_traits), 
                                              SM_QTL = SM_e,
                                              SM_Observed = SM,
                                              LG_QTL = LG_e,
                                              F3_Observed = F3,
                                              LG_Observed = LG)) %>% separate(variable, c("Line", "Type"))
growth_pred_plot_SUR = ggplot() + scale_x_discrete(labels = paste("Week", 1:7)) + labs(y = "Weekly growth (g)", x = "Start week") + geom_line(data = post, color = "gray", size = 0.5, linetype = 1, alpha = 0.1, aes(trait, value, group = interaction(Type, Line, iterations))) + geom_line(size = 1, data = growth_prediction, aes(trait, value, group = interaction(Type, Line), color = Line, linetype = Type)) + theme_cowplot()
save_plot("data/growth_multiple_regression_ancestral_prediction.png", growth_pred_plot_SUR, base_height = 7, base_aspect_ratio = 2)

growth_observed_plot = ggplot() + scale_x_discrete(labels = paste("Week", 1:7)) + labs(y = "Weekly growth (g)", x = "Start week") + geom_line(size = 1, data = filter(growth_prediction, Type == "Observed"), aes(trait, value, group = Line, color = Line)) + theme_cowplot()
save_plot("data/growth_LG_SM_F3.png", growth_observed_plot, base_height = 7, base_aspect_ratio = 2)


par(mar = c(0, 0, 0, 0))
g_plot = ~corrplot.mixed(cov2cor(G),       upper = "ellipse", mar = c(0, 0, 0, 0))
plot_grid(growth_observed_plot, g_plot, labels = LETTERS[1:2])

growth_observed_parentals_plot = ggplot() + scale_x_discrete(labels = paste("Week", 1:7)) + labs(y = "Weekly growth (g)", x = "Start week") + geom_line(size = 1, data = filter(growth_prediction, Type == "Observed", Line != "F3"), aes(trait, value, group = Line, color = Line)) + scale_color_manual(values=c("#00BA38", "#619CFF")) + theme_cowplot()
save_plot("data/growth_LG_SM.png", growth_observed_parentals_plot, base_height = 7, base_aspect_ratio = 2)

growth_observed_parentals_Dz_plot = ggplot() + scale_x_discrete(labels = paste("Week", 1:7)) + labs(y = "Weekly growth (g)", x = "Start week") + geom_line(size = 1, data = filter(growth_prediction, Type == "Observed", Line != "F3"), aes(trait, value, group = Line, color = Line)) + scale_color_manual(values=c("#00BA38", "#619CFF")) + 
  annotate("segment", x = 1:7, xend = 1:7, y = SM, yend = LG) + 
  annotate("text", x = 3.5, y = 6, label = "Phenotypic\ndivergence") + theme_cowplot()
save_plot("data/growth_LG_SM_DZ1.png", growth_observed_parentals_Dz_plot, base_height = 7, base_aspect_ratio = 2)

growth_observed_parentals_Dz2_plot = ggplot() + scale_x_discrete(labels = paste("Week", 1:7)) + labs(y = "Weekly growth (g)", x = "Start week") + geom_line(size = 1, data = filter(growth_prediction, Type == "Observed", Line != "F3"), aes(trait, value, group = Line, color = Line)) + scale_color_manual(values=c("#00BA38", "#619CFF")) +
  geom_line(data = data.frame(trait = as.factor(growth_traits), dz = d_z), aes(trait, dz, group = 1)) +
  annotate("text", x = 3.5, y = 5.5, label = "Phenotypic\ndivergence")
save_plot("data/growth_LG_SM_DZ2.png", growth_observed_parentals_Dz2_plot, base_height = 7, base_aspect_ratio = 2)

library(scales)
show_col(hue_pal()(3))


w_ad = rstan::extract(full_HCp, pars = "w_ad")[[1]]
effect_matrix = aaply(w_ad, c(2, 3), mean)
SM_e = rowSums(-1 * as.matrix(effect_matrix)) * sapply(growth_phen[,growth_traits], sd) + 
    sapply(growth_phen[,growth_traits], mean)
LG_e = rowSums(as.matrix(effect_matrix)) * sapply(growth_phen[,growth_traits], sd) + 
    sapply(growth_phen[,growth_traits], mean)
post = adply(w_ad, 1, posterior_predict)
growth_prediction = reshape2::melt(data.frame(trait = as.factor(growth_traits), 
                                              SM_GenPred = SM_e,
                                              SM_Observed = SM,
                                              LG_GenPred = LG_e,
                                              F3_Observed = F3,
                                              LG_Observed = LG)) %>% separate(variable, c("Line", "Type"))
growth_pred_plot_HC_full = ggplot() + scale_x_discrete(labels = paste("Week", 1:7)) + labs(y = "Weekly growth (g)", x = "Start week") + geom_line(data = post, color = "gray", size = 0.5, linetype = 1, alpha = 0.1, aes(trait, value, group = interaction(Type, Line, iterations))) + geom_line(size = 1, data = growth_prediction, aes(trait, value, group = interaction(Type, Line), color = Line, linetype = Type))
save_plot("data/growth_multiple_regression_ancestral_prediction_full_genome.png", growth_pred_plot_HC_full, base_height = 7, base_aspect_ratio = 2)

plot_grid(growth_pred_plot_SUR, growth_pred_plot_HC_full)


## Pleiotropic partition

pleiotropic_partition = a_effect_matrix
pleiotropic_partition[growth_traits] = (pleiotropic_partition[growth_traits]^2)
pleiotropic_partition$id = factor(pleiotropic_partition$id, levels = pleiotropic_partition$id)
a_pleiotropic_partition = pleiotropic_partition %>% 
  separate(id, c("chrom", "marker")) %>% 
  mutate(chrom = as.numeric(chrom), 
         marker = as.numeric(marker)) %>%
  arrange(chrom, marker) %>%
  mutate(id = factor(paste(chrom, marker, sep= "_"), levels = paste(chrom, marker, sep= "_")))

a_pleiotropic_partition_plot =  
  ggplot() + 
  geom_bar(data = gather(a_pleiotropic_partition, key, value, growth_traits), 
           aes(id, value, group = key, color = key, fill = key), stat = "identity") +
  scale_y_continuous(limits = c(0, 0.25)) + background_grid(major = "xy", minor = "none") + 
  scale_fill_viridis(discrete = TRUE, option = "D", 
                     guide = guide_legend(direction = "horizontal", 
                                          label.position = "top",
                                          nrow = 1, title = NULL)) + 
  scale_color_viridis(discrete = TRUE, option = "D", 
                      guide = guide_legend(direction = "horizontal",
                                           title = NULL)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = c(0.5, 0.9)) + 
  labs(x = "Marker", y = "Squared contribution to\n scaled pleitropic vector")

pleiotropic_partition = d_effect_matrix 
pleiotropic_partition[growth_traits] = (pleiotropic_partition[growth_traits]^2)
pleiotropic_partition$id = factor(pleiotropic_partition$id, levels = pleiotropic_partition$id)
d_pleiotropic_partition = pleiotropic_partition %>% 
  separate(id, c("chrom", "marker")) %>% 
  mutate(chrom = as.numeric(chrom), 
         marker = as.numeric(marker)) %>%
  arrange(chrom, marker) %>%
  mutate(id = factor(paste(chrom, marker, sep= "_"), levels = paste(chrom, marker, sep= "_")))
d_pleiotropic_partition_plot =  
  ggplot() + 
  geom_bar(data = gather(d_pleiotropic_partition, key, value, growth_traits), aes(id, value, group = key, color = key, fill = key), stat = "identity") +
  scale_fill_viridis(discrete = TRUE, option = "D") + scale_color_viridis(discrete = TRUE, option = "D") + 
  scale_y_continuous(limits = c(0, 0.25)) + background_grid(major = "xy", minor = "none") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  labs(x = "Marker", y = "Squared contribution to\n scaled pleitropic vector")
pleiotropic_Effects_ad_dm = plot_grid(a_pleiotropic_partition_plot, d_pleiotropic_partition_plot, ad_biplot, ncol = 1, labels = c("A", "B"))
save_plot("data/growth_pleiotropic_partition_ad_dm.png", pleiotropic_Effects_ad_dm, base_height = 4.5, base_aspect_ratio = 2.5, nrow = 3) 


## Pleiotropic partition Genome prediction

pleiotropic_partition = a_effect_matrix_HC
pleiotropic_partition[growth_traits] = sqrt(pleiotropic_partition[growth_traits]^2)
pleiotropic_partition$id = factor(pleiotropic_partition$id, levels = pleiotropic_partition$id)
a_pleiotropic_partition = pleiotropic_partition %>% 
  separate(id, c("chrom", "marker")) %>% 
  mutate(chrom = as.numeric(chrom), 
         marker = as.numeric(marker)) %>%
  arrange(chrom, marker) %>%
  mutate(id = factor(paste(chrom, marker, sep= "_"), levels = paste(chrom, marker, sep= "_")))

a_pleiotropic_partition_plot =  
  ggplot() + 
  geom_bar(data = gather(a_pleiotropic_partition, key, value, growth_traits), 
           aes(id, value, group = key, color = key, fill = key), stat = "identity") +
  scale_y_continuous(limits = c(0, 0.35)) + background_grid(major = "xy", minor = "none") + 
  scale_fill_viridis(discrete = TRUE, option = "D", 
                     guide = guide_legend(direction = "horizontal", 
                                          label.position = "top",
                                          nrow = 1, title = NULL)) + 
  scale_color_viridis(discrete = TRUE, option = "D", 
                      guide = guide_legend(direction = "horizontal",
                                           title = NULL)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = c(0.45, 0.9)) + 
  labs(x = "Marker", y = "Squared contribution to\n scaled pleitropic vector")
#save_plot("data/growth_pleiotropic_partition_additive.png", pleiotropic_partition_plot, base_height = 7, base_aspect_ratio = 2) 

pleiotropic_partition = d_effect_matrix_HC 
pleiotropic_partition[growth_traits] = sqrt(pleiotropic_partition[growth_traits]^2)
pleiotropic_partition$id = factor(pleiotropic_partition$id, levels = pleiotropic_partition$id)
d_pleiotropic_partition = pleiotropic_partition %>% 
  separate(id, c("chrom", "marker")) %>% 
  mutate(chrom = as.numeric(chrom), 
         marker = as.numeric(marker)) %>%
  arrange(chrom, marker) %>%
  mutate(id = factor(paste(chrom, marker, sep= "_"), levels = paste(chrom, marker, sep= "_")))
d_pleiotropic_partition_plot =  
  ggplot() + 
  geom_bar(data = gather(d_pleiotropic_partition, key, value, growth_traits), aes(id, value, group = key, color = key, fill = key), stat = "identity") +
  scale_fill_viridis(discrete = TRUE, option = "D") + scale_color_viridis(discrete = TRUE, option = "D") + 
  scale_y_continuous(limits = c(0, 0.35)) + background_grid(major = "xy", minor = "none") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  labs(x = "Marker", y = "Squared contribution to\n scaled pleitropic vector")
pleiotropic_Effects_ad_dm = plot_grid(a_pleiotropic_partition_plot, d_pleiotropic_partition_plot, ad_biplot_HC, ncol = 1, labels = c("A", "B"))
save_plot("data/growth_pleiotropic_partition_ad_dm_GP.png", pleiotropic_Effects_ad_dm, base_height = 4.5, base_aspect_ratio = 2.5, nrow = 3) 
#save_plot("data/growth_pleiotropic_partition_dominance.png", pleiotropic_partition_plot, base_height = 7, base_aspect_ratio = 2) 