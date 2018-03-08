if(!require(purrr)){install.packages("purrr"); library(purrr)}
if(!require(corrplot)){install.packages("corrplot"); library(corrplot)}

setwd("/home/diogro/projects/mouse-qtls")
source('read_mouse_data.R')

#Rdatas_folder = "~/gdrive/LGSM_project_Rdatas/"
Rdatas_folder = "./Rdatas/"

vectorCor = function(x, y) Normalize(x) %*% Normalize(y)
markerCov = function(marker1, marker2){
  marker1_col = growth_markers[,makeMarkerList(marker1)]
  marker2_col = growth_markers[,makeMarkerList(marker2)]
  cov(marker1_col, marker2_col)[1]
}
makeMarkerList = function(pos) paste('chrom', pos[1],"_", 'A', pos[2], sep = '')
lt = function(x) x[lower.tri(x, diag = TRUE)]

install_load("doMC", "lme4qtl", "qvalue")
registerDoMC(15)

growth_data = inner_join(growth_phen_std,
                        growth_markers,
                        by = "ID") %>%
  gather(variable, value, growth12:growth78)
growth_data = mutate(growth_data, FAMILY = as.factor(FAMILY), variable = as.factor(variable))

e = read_csv("./data/growth_significant_marker_effects_SUR.csv") %>% arrange(trait, id, class)
load(file = "./Rdatas/significant_stan_fit.Rdata")
eHC = read_csv("./data/growth_significant_marker_effectsHC.csv") %>% arrange(trait, id, class)
load(file = "./Rdatas/significant_stan_HCp_fit.Rdata")
full_HCp = readRDS("./Rdatas/growth_scaled_allmarkers_HCPlus")

a_effect_matrix = e %>%
    select(id, class, trait, mean) %>%
    spread(trait, mean) %>%
    filter(class == "additive") %>% select(-class)
a_effect_matrix_HC = eHC %>%
    select(id, class, trait, mean) %>%
    spread(trait, mean) %>%
    filter(class == "additive") %>% select(-class)
d_effect_matrix = e %>%
  select(id, class, trait, mean) %>%
  spread(trait, mean) %>%
  filter(class == "dominance") %>% select(-class)
d_effect_matrix_HC = eHC %>%
  select(id, class, trait, mean) %>%
  spread(trait, mean) %>%
  filter(class == "dominance") %>% select(-class)

png("data/growth_additive_effects_PCA.png")
biplot(prcomp(a_effect_matrix[,growth_traits]))
dev.off()
png("data/growth_dominance_effects_PCA.png")
biplot(prcomp(d_effect_matrix[,growth_traits]))
dev.off()
eig_effects  =  eigen(cov(effect_matrix[,growth_traits]))
plot(eig_effects$values)
abline(0, 1)


LG = c(3.785,4.435,8.43,7.395,2.995,1.85,2.085)
SM = c(3.31 ,2.98,3.82,2.175,0.765,1.165,0.51)
F3 = sapply(growth_phen[,growth_traits], mean)

d_z = LG - SM
library(evolqg)
load("./Rdatas/growth_CovMatrices.Rdata")
growth_sds = apply(growth_phen[,growth_traits], 2, sd)
G = G_stan #* outer(growth_sds, growth_sds)
corrplot.mixed(cov2cor(G), upper = "ellipse")
plot(eigen(G_mcmc)$values)
G_ext4 = ExtendMatrix(G, ret.dim = 4)[[1]]
G_ext5 = ExtendMatrix(G, ret.dim = 5)[[1]]
solve(G, d_z)
(beta = solve(G_ext4, d_z))
solve(G_ext5, d_z)


mean_a = colMeans(a_effect_matrix[,growth_traits])
mean_d = colMeans(d_effect_matrix[,growth_traits])
vectorCor(mean_a, d_z)
vectorCor(mean_d, d_z)



markerMatrix
calcVa = function(i, a_effects, markerMatrix){
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
  # Variance due to LD with focal marker
  for(j in 1:nrow(markerMatrix)){
    if (i != j) V_a = V_a + markerCov(focal_marker, markerMatrix[j,]) * 
                            outer(effects[,i], effects[,j])
  }
  V_a
}
calcVd = function(i, d_effects, markerMatrix){
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
calcVa(1, effect_matrix, significantMarkerMatrix)

w_ad = rstan::extract(stan_model_SUR, pars = "w_ad")[[1]]
effect_matrix_additive = aaply(w_ad, c(2, 3), mean)
w_dm = rstan::extract(stan_model_SUR, pars = "w_dm")[[1]]
effect_matrix_dominance = aaply(w_dm, c(2, 3), mean)
save(w_ad, w_dm, effect_matrix_additive, effect_matrix_dominance, file = "Rdatas/growth_add_dom_effectsMatrix.Rdata")
  
Va = laply(seq_along(1:400), function(i) colSums(laply(seq_along(significantMarkerList), calcVa, w_ad[i,,], significantMarkerMatrix)), .parallel = TRUE)

Vd = laply(seq_along(1:400), function(i) colSums(laply(seq_along(significantMarkerList), calcVd, w_dm[i,,], significantMarkerMatrix)), .parallel = TRUE)

Va_mean = aaply(Va, c(2, 3), mean)
Va_upper = aaply(Va, c(2, 3), quantile, 0.975)
Va_lower = aaply(Va, c(2, 3), quantile, 0.025)

Vd_mean  = aaply(Vd, c(2, 3), mean)
Vd_lower = aaply(Vd, c(2, 3), quantile, 0.025)
Vd_upper = aaply(Vd, c(2, 3), quantile, 0.975)

G = aaply(Gs_stan, c(2, 3), mean)
G_lower = aaply(Gs_stan, c(2, 3), quantile, 0.025)
G_upper = aaply(Gs_stan, c(2, 3), quantile, 0.975)

Vg = 1/2 * Va + 1/4 * Vd

Vg_mean  = aaply(Vg, c(2, 3), mean)
Vg_lower = aaply(Vg, c(2, 3), quantile, 0.025)
Vg_upper = aaply(Vg, c(2, 3), quantile, 0.975)

corrplot.mixed(G, upper = "ellipse")
corrplot.mixed(Vg_mean, upper = "ellipse")
write.csv(MatrixCompare(Vd_mean, G), file = "./data/Vd_FamilyG_comparison.csv")
write.csv(MatrixCompare(Va_mean, G), file = "./data/Va_FamilyG_comparison.csv")
write.csv(MatrixCompare(Vg_mean, G), file = "./data/Vg_FamilyG_comparison.csv")

data.frame(Vg = diag(Vg_mean), G = diag(G)) %>% gather %>%
  ggplot(aes(c(1:7, 1:7), value, group = key, color = key)) + geom_point() + geom_line()
png("./data/growth_family_qtl_cov.png", width = 1080, height = 640)
plot(lt(G)~lt(Vg_mean), pch = 19, 
     ylab = "Family G-matrix covariances", xlab = "Genetic covariances predicted from QTLs (1/2 * Va + 1/4 * Vd)", 
     main = "Growth traits", xlim = c(-0.03, 0.17), ylim = c(-0.22, 0.6))
segments(x0 = lt(Vg_lower), y0 = lt(G), x1 = lt(Vg_upper), y1 = lt(G))
segments(x0 = lt(Vg_mean), y0 = lt(G_lower), x1 = lt(Vg_mean), y1 = lt(G_upper))
points(diag(G)~diag(Vg_mean), col = "tomato3", pch = 19)
abline(lm(lt(G)~lt(Vg_mean)))
abline(0, 1, col = "blue")
text(0.15, 0.12, "Identity", col = "blue")
text(0.11, 0.35, "Variances", col = "tomato3")
text(0.05, -0.05, "Co-variances")
dev.off()
summary(lm(lt(G)~lt(Vg_mean)))

w_ad = rstan::extract(full_HCp, pars = "w_ad")[[1]]
effect_matrix_HC = aaply(w_ad, c(2, 3), mean)
Va_HC = colSums(laply(seq_along(significantMarkerList), calcVa, effect_matrix_HC, markerMatrix))

library(viridis)

vectorCor(d_z, beta)
random_vec = matrix(rnorm(7*1000), 1000, 7)
quantile(abs(apply(random_vec, 1, vectorCor, rep(1, 7))), 0.95)
crss = data.frame(beta = apply(a_effect_matrix[,growth_traits], 1, vectorCor, beta),
                    dz = abs(apply(a_effect_matrix[,growth_traits], 1, vectorCor, d_z))) %>% gather(class, value, beta:dz)

corrs = data.frame(beta = apply(a_effect_matrix[,growth_traits], 1, vectorCor, beta),
                     dz = abs(apply(a_effect_matrix[,growth_traits], 1, vectorCor, d_z)),
                   norm = apply(a_effect_matrix[,growth_traits], 1, Norm))
ggplot(crss, aes(class, value, fill = class)) + geom_violin()
ggplot(corrs, aes(norm, beta)) + geom_point() + geom_smooth(method = "lm", color = "black", se = FALSE)
ggplot(corrs, aes(dz, norm)) + geom_point() + geom_smooth(method = "lm", color = "black", se = FALSE)

lm(beta~norm, data = corrs) %>% summary
lm(norm~dz, data = corrs) %>% summary

growth_m = as.numeric(ddply(growth_phen, .(SEX), numcolwise(mean))[2,growth_traits])
growth_f = as.numeric(ddply(growth_phen, .(SEX), numcolwise(mean))[1,growth_traits])

growth_prediction = reshape2::melt(data.frame(trait = as.factor(growth_traits), 
                                              SM_SUR = SM_e,
                                              SM_HC = SM_e_HC,
                                              SM_Observed = SM,
                                              LG_SUR = LG_e,
                                              LG_HC = LG_e_HC,
                                              F3_Observed = F3,
                                              LG_Observed = LG)) %>% separate(variable, c("Line", "Type"))
growth_pred_plot = ggplot(growth_prediction, aes(trait, value, group = interaction(Type, Line), color = Line, linetype = Type)) + geom_line(size = 1) + scale_x_discrete(labels = paste("Week", 1:7)) + labs(y = "Weekly growth (g)", x = "Start week")


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
                                              SM_stan = SM_e,
                                              SM_Observed = SM,
                                              LG_stan = LG_e,
                                              F3_Observed = F3,
                                              LG_Observed = LG)) %>% separate(variable, c("Line", "Type"))
growth_pred_plot_SUR = ggplot() + scale_x_discrete(labels = paste("Week", 1:7)) + labs(y = "Weekly growth (g)", x = "Start week") + geom_line(data = post, color = "gray", size = 0.5, linetype = 1, alpha = 0.1, aes(trait, value, group = interaction(Type, Line, iterations))) + geom_line(size = 1, data = growth_prediction, aes(trait, value, group = interaction(Type, Line), color = Line, linetype = Type))
save_plot("data/growth_multiple_regression_ancestral_prediction.png", growth_pred_plot_SUR, base_height = 7, base_aspect_ratio = 2)

w_ad = rstan::extract(full_HCp, pars = "w_ad")[[1]]
effect_matrix = aaply(w_ad, c(2, 3), mean)
SM_e = rowSums(-1 * as.matrix(effect_matrix)) * sapply(growth_phen[,growth_traits], sd) + 
    sapply(growth_phen[,growth_traits], mean)
LG_e = rowSums(as.matrix(effect_matrix)) * sapply(growth_phen[,growth_traits], sd) + 
    sapply(growth_phen[,growth_traits], mean)
post = adply(w_ad, 1, posterior_predict)
growth_prediction = reshape2::melt(data.frame(trait = as.factor(growth_traits), 
                                              SM_stan = SM_e,
                                              SM_Observed = SM,
                                              LG_stan = LG_e,
                                              F3_Observed = F3,
                                              LG_Observed = LG)) %>% separate(variable, c("Line", "Type"))
growth_pred_plot_HC_full = ggplot() + scale_x_discrete(labels = paste("Week", 1:7)) + labs(y = "Weekly growth (g)", x = "Start week") + geom_line(data = post, color = "gray", size = 0.5, linetype = 1, alpha = 0.1, aes(trait, value, group = interaction(Type, Line, iterations))) + geom_line(size = 1, data = growth_prediction, aes(trait, value, group = interaction(Type, Line), color = Line, linetype = Type))
plot_grid(growth_pred_plot_SUR, growth_pred_plot_HC_full)

effect_matrix[,1]
Reduce("+", alply(effect_matrix, 2, function(x) outer(x, x)))
