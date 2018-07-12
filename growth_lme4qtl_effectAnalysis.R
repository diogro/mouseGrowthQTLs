#setwd("/home/diogro/projects/mouse-qtls")
source('read_mouse_data.R')

ncores = 8
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
G = G_stan
#png("~/G_LGSM", width = 600, height = 600)
corrplot.mixed(cov2cor(G), upper = "ellipse")
#dev.off()
plot(eigen(G)$values)
G_ext4 = ExtendMatrix(G, ret.dim = 4)[[1]]
G_ext5 = ExtendMatrix(G, ret.dim = 5)[[1]]
solve(G, d_z)
(beta = Normalize(solve(G_ext4, d_z)))
beta_w = Normalize((G_stan_w %*% c(1, rep(0, length(growth_traits))))[-1])
solve(G_ext5, d_z)

scale = Norm(G %*% beta_w)/Norm(d_z)
beta_scaled = (beta_w/scale)
d_z_w = (G %*% beta_scaled)[,1]
vectorCor(d_z_w, d_z)
Norm(d_z_w)
Norm(d_z)

mean_a = colMeans(aaply(as.matrix(a_effect_matrix[,growth_traits]), 1, function(x) x * growth_phen_sd))
mean_d = colMeans(aaply(as.matrix(d_effect_matrix[,growth_traits]), 1, function(x) x * growth_phen_sd))
vectorCor(mean_a, d_z)
vectorCor(mean_d, d_z)

vectorCor(mean_a, beta_w)
vectorCor(mean_d, beta_w)
vectorCor(beta_w, d_z)
vectorCor(beta, d_z)
vectorCor(beta_w, beta)

calcVa = function(i, a_effects, d_effects, markerMatrix, include_LD = TRUE){
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
  if(include_LD){
    for(j in 1:nrow(markerMatrix)){
      if (i != j) V_a = V_a + markerCov(focal_marker, markerMatrix[j,]) * 
          outer(a_effects[,i], a_effects[,j])
    }
  }
  V_a
}
calcVd = function(i, d_effects, markerMatrix, include_LD = TRUE){
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
  if(include_LD){
    for(j in 1:nrow(markerMatrix)){
      d2ij = markerCov(focal_marker, markerMatrix[j,])
      if (i != j) V_d = V_d + d2ij^2 * outer(d_effects[,i], d_effects[,j])
    }
  }
  V_d
}

w_ad = rstan::extract(stan_model_SUR, pars = "w_ad")[[1]]
effect_matrix_additive = aaply(w_ad, c(2, 3), mean)
w_dm = rstan::extract(stan_model_SUR, pars = "w_dm")[[1]]
effect_matrix_dominance = aaply(w_dm, c(2, 3), mean)
save(w_ad, w_dm, effect_matrix_additive, effect_matrix_dominance, file = "Rdatas/growth_add_dom_effectsMatrix.Rdata")
load(file = "Rdatas/growth_add_dom_effectsMatrix.Rdata")

 # Va = laply(seq_along(1:400), function(i) colSums(laply(1:nrow(significantMarkerMatrix), calcVa, w_ad[i,,], w_dm[i,,], significantMarkerMatrix)), .parallel = TRUE)
 # Va_pleio = laply(seq_along(1:400), function(i) colSums(laply(1:nrow(significantMarkerMatrix), calcVa, w_ad[i,,], w_dm[i,,], significantMarkerMatrix, include_LD = FALSE)), .parallel = TRUE)
 # Vd = laply(seq_along(1:400), function(i) colSums(laply(1:nrow(significantMarkerMatrix), calcVd, w_dm[i,,], significantMarkerMatrix)), .parallel = TRUE)
 # Vd_pleio = laply(seq_along(1:400), function(i) colSums(laply(1:nrow(significantMarkerMatrix), calcVd, w_dm[i,,], significantMarkerMatrix, include_LD = FALSE)), .parallel = TRUE)
 # save(Va,Va_pleio, Vd, Vd_pleio, file = paste0(Rdatas_folder, "VaVd_QTL.Rdata"))
load(paste0(Rdatas_folder, "VaVd_QTL.Rdata"))

Va_mean = aaply(Va, c(2, 3), mean)
Va_upper = aaply(Va, c(2, 3), quantile, 0.975)
Va_lower = aaply(Va, c(2, 3), quantile, 0.025)

VaP_mean = aaply(Va_pleio, c(2, 3), mean)
VaP_upper = aaply(Va_pleio, c(2, 3), quantile, 0.975)
VaP_lower = aaply(Va_pleio, c(2, 3), quantile, 0.025)

Vd_mean  = aaply(Vd, c(2, 3), mean)
Vd_lower = aaply(Vd, c(2, 3), quantile, 0.025)
Vd_upper = aaply(Vd, c(2, 3), quantile, 0.975)

VdP_mean  = aaply(Vd_pleio, c(2, 3), mean)
VdP_lower = aaply(Vd_pleio, c(2, 3), quantile, 0.025)
VdP_upper = aaply(Vd_pleio, c(2, 3), quantile, 0.975)

G = aaply(Gs_stan, c(2, 3), mean)
dimnames(G) = list(1:7, 1:7)
G_lower = aaply(Gs_stan, c(2, 3), quantile, 0.025)
G_upper = aaply(Gs_stan, c(2, 3), quantile, 0.975)

#G_dam_lower = aaply(Gs_dam, c(2, 3), quantile, 0.025)
#G_dam_upper = aaply(Gs_dam, c(2, 3), quantile, 0.975)

Vg = 0.5 * Va + 0.25 * Vd

Vg_mean  = aaply(Vg, c(2, 3), mean)
Vg_lower = aaply(Vg, c(2, 3), quantile, 0.025)
Vg_upper = aaply(Vg, c(2, 3), quantile, 0.975)

Vg_pleio = 0.5 * Va_pleio + 0.25 * Vd_pleio

VgP_mean  = aaply(Vg_pleio, c(2, 3), mean)
VgP_lower = aaply(Vg_pleio, c(2, 3), quantile, 0.025)
VgP_upper = aaply(Vg_pleio, c(2, 3), quantile, 0.975)

matrices <- list(FullSib = G_stan,
                 "Va QTL" = Va_mean,
                 "Vd QTL" = Vd_mean,
                 "1/2 Va + 1/4 Vd" = Vg_mean,
                 "FullSib cross-foster" = G_cf,
                 "FullSib non-cross-foster" = G_ncf,
                 #"G_dam" = G_dam,
                 #"G_nurse" = G_nurse,
                 "beta_r" = beta,
                 "beta_w" = beta_w,
                 "delta Z" = d_z)
PrintMatrix(matrices, output.file = "./matrix.csv")

png("./data/growth_family_Vg_FullSibG_correlation.png", width = 1500, height = 1500)
old.par = par()
par(mfrow = c(2, 2), cex=2, oma = c(0, 0, 0, 0))
corrplot.mixed(cov2cor(G),       upper = "ellipse", mar = c(0, 0, 1, 0))
corrplot.mixed(cov2cor(Vg_mean), upper = "ellipse", mar = c(0, 0, 1, 0))
corrplot.mixed(cov2cor(Va_mean), upper = "ellipse", mar = c(0, 0, 1, 0))
corrplot.mixed(cov2cor(Vd_mean), upper = "ellipse", mar = c(0, 0, 1, 0))
dev.off()
par(old.par)

png("./data/growth_Va_PleiotropyLD_and_OnlyPleitropy.png", width = 1500, height = 800)
par(mfrow = c(1, 2), cex=2, oma = c(0, 0, 0, 0))
corrplot.mixed(cov2cor(Va_mean), upper = "ellipse", mar = c(0, 0, 1, 0), main = "Pleiotropy + LD")
corrplot.mixed(cov2cor(VaP_mean), upper = "ellipse", mar = c(0, 0, 1, 0), main = "Pleiotropy only")
dev.off()

diag(VaP_mean)
diag(Va_mean)
MatrixCompare(Va_mean, VaP_mean)

write.csv(MatrixCompare(Vd_mean, G), file = "./data/TalkStuff/Vd_FamilyG_comparison.csv")
write.csv(MatrixCompare(Va_mean, G), file = "./data/TalkStuff/Va_FamilyG_comparison.csv")
write.csv(MatrixCompare(Vg_mean, G), file = "./data/TalkStuff/Vg_FamilyG_comparison.csv")
write.csv(MatrixCompare(Vg_mean, Vd_mean), file = "./data/TalkStuff/Va_Vd_comparison.csv")

# write.csv(MatrixCompare(Vd_mean, G_dam), file = "./data/TalkStuff/Vd_DamG_comparison.csv")
# write.csv(MatrixCompare(Va_mean, G_dam), file = "./data/TalkStuff/Va_DamG_comparison.csv")
# write.csv(MatrixCompare(Vg_mean, G_dam), file = "./data/TalkStuff/Vg_DamG_comparison.csv")
# write.csv(MatrixCompare(G, G_dam), file = "./data/TalkStuff/GFamily_GDam_comparison.csv")

# write.csv(MatrixCompare(Vd_mean, G_nurse), file = "./data/TalkStuff/Vd_nurseG_comparison.csv")
# write.csv(MatrixCompare(Va_mean, G_nurse), file = "./data/TalkStuff/Va_nurseG_comparison.csv")
# write.csv(MatrixCompare(Vg_mean, G_nurse), file = "./data/TalkStuff/Vg_nurseG_comparison.csv")
# write.csv(MatrixCompare(G, G_nurse), file = "./data/TalkStuff/GFamily_Gnurse_comparison.csv")
#write.csv(MatrixCompare(G_dam, G_nurse), file = "./data/TalkStuff/Gdam_Gnurse_comparison.csv")

g_predict = function() {plot(lt(G)~lt(Vg_mean), pch = 19, 
     ylab = "Mixed model\nG-matrix", xlab = "Genetic covariances predicted\n from pleiotropic vectors", 
     main = "", xlim = c(-0.22, 0.55), ylim = c(-0.22, 0.55))
segments(x0 = lt(Vg_lower), y0 = lt(G), x1 = lt(Vg_upper), y1 = lt(G))
segments(x0 = lt(Vg_mean), y0 = lt(G_lower), x1 = lt(Vg_mean), y1 = lt(G_upper))
points(diag(G)~diag(Vg_mean), col = "tomato3", pch = 19)
abline(lm(lt(G)~lt(Vg_mean)))
abline(0, 1, col = "blue")
abline(v = 0)
abline(h = 0)
text(0.2, 0.12, "Identity", col = "blue")
text(0.25, 0.35, "Variances", col = "tomato3")
text(0.1, -0.05, "Co-variances")}; g_predict()

summary(lm(lt(G)~lt(Vg_mean)))

line = 2.5
cex = 2
las =2
padj = -10.5
png("data/growth_cov_prediction_composite.png", width = 900, height = 900)
par(cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
layout(matrix(c(1,2,3, 4), 2, 2, byrow = TRUE))
corrplot.mixed(cov2cor(Va_mean), upper = "ellipse", mar = c(0, 0, 2, 0), main = "Additive correlations")
mtext("A", side=2, line=line, cex=cex, las=las, padj = padj)
corrplot.mixed(cov2cor(Vd_mean), upper = "ellipse", mar = c(0, 0, 2, 0), main = "Dominance correlations")
mtext("B", side=2, line=line, cex=cex, las=las, padj = padj)
corrplot.mixed(cov2cor(Vg_mean), upper = "ellipse", mar = c(0, 0, 2, 0), main = "Genetic correlations")
mtext("C", side=2, line=line, cex=cex, las=las, padj = padj)
par(mar = c(6, 7, 1, 1), mgp=c(4,1,0))
g_predict()
mtext("D", side=2, line=line, cex=cex, las=las, padj = padj)
dev.off()
# png("./data/growth_dam_qtl_cov.png", width = 1500, height = 800)
# par(mfrow = c(1, 1), cex=2)
# plot(lt(G_dam)~lt(Vg_mean), pch = 19, 
#      ylab = "Dam G-matrix covariances", xlab = "Genetic covariances predicted from QTLs (1/2 * Va + 1/4 * Vd)", 
#      main = "Growth traits", xlim = c(-0.03, 0.17), ylim = c(-0.22, 0.6))
# segments(x0 = lt(Vg_lower), y0 = lt(G_dam), x1 = lt(Vg_upper), y1 = lt(G_dam))
# segments(x0 = lt(Vg_mean), y0 = lt(G_dam_lower), x1 = lt(Vg_mean), y1 = lt(G_dam_upper))
# points(diag(G_dam)~diag(Vg_mean), col = "tomato3", pch = 19)
# abline(lm(lt(G_dam)~lt(Vg_mean)))
# abline(0, 1, col = "blue")
# text(0.15, 0.12, "Identity", col = "blue")
# text(0.11, 0.35, "Variances", col = "tomato3")
# text(0.05, -0.05, "Co-variances")
# dev.off()

# Va_GP = laply(seq_along(1:400), function(i) colSums(laply(1:nrow(significantMarkerMatrix),
#                                                           calcVa,
#                                                           w_ad_HC[i,,],
#                                                           w_dm_HC[i,,],
#                                                           markerMatrix)),
#               .parallel = TRUE)
# Va_GP_pleio = laply(seq_along(1:400), function(i) colSums(laply(1:nrow(significantMarkerMatrix),
#                                                           calcVa,
#                                                           w_ad_HC[i,,],
#                                                           w_dm_HC[i,,],
#                                                           markerMatrix, include_LD = FALSE)),
#               .parallel = TRUE)
# Vd_GP = laply(seq_along(1:400), function(i) colSums(laply(1:nrow(significantMarkerMatrix),
#                                                           calcVd,
#                                                           w_dm_HC[i,,],
#                                                           markerMatrix)),
#               .parallel = TRUE)
# Vd_GP_pleio = laply(seq_along(1:400), function(i) colSums(laply(1:nrow(significantMarkerMatrix),
#                                                           calcVd,
#                                                           w_dm_HC[i,,],
#                                                           markerMatrix, include_LD = FALSE)),
#               .parallel = TRUE)
# 
# save(Va_GP, Va_GP_pleio, Vd_GP, Vd_GP_pleio, file = paste0(Rdatas_folder, "VaVd_GP_QTL.Rdata"))
load(paste0(Rdatas_folder, "VaVd_GP_QTL.Rdata"))

Va_GP_mean = aaply(Va_GP, c(2, 3), mean)
Va_GP_upper = aaply(Va_GP, c(2, 3), quantile, 0.975)
Va_GP_lower = aaply(Va_GP, c(2, 3), quantile, 0.025)

Va_GP_pleio_mean = aaply(Va_GP_pleio, c(2, 3), mean)
Va_GP_pleio_upper = aaply(Va_GP_pleio, c(2, 3), quantile, 0.975)
Va_GP_pleio_lower = aaply(Va_GP_pleio, c(2, 3), quantile, 0.025)

Vd_GP_mean  = aaply(Vd_GP, c(2, 3), mean)
Vd_GP_lower = aaply(Vd_GP, c(2, 3), quantile, 0.025)
Vd_GP_upper = aaply(Vd_GP, c(2, 3), quantile, 0.975)

Vd_GP_pleio_mean  = aaply(Vd_GP_pleio, c(2, 3), mean)
Vd_GP_pleio_lower = aaply(Vd_GP_pleio, c(2, 3), quantile, 0.025)
Vd_GP_pleio_upper = aaply(Vd_GP_pleio, c(2, 3), quantile, 0.975)

Vg_GP = 0.5 * Va_GP + 0.25 * Vd_GP

Vg_GP_mean  = aaply(Vg_GP, c(2, 3), mean)
Vg_GP_lower = aaply(Vg_GP, c(2, 3), quantile, 0.025)
Vg_GP_upper = aaply(Vg_GP, c(2, 3), quantile, 0.975)

g_predict_GP = function() {
  plot(lt(G)~lt(Vg_GP_mean), pch = 19, 
                             ylab = "Mixed model\nG-matrix", xlab = "Genetic covariances predicted\n from pleiotropic vectors", 
                             main = "", xlim = c(-0.22, 0.55), ylim = c(-0.22, 0.55))
  segments(x0 = lt(Vg_GP_lower), y0 = lt(G), x1 = lt(Vg_GP_upper), y1 = lt(G))
  segments(x0 = lt(Vg_GP_mean), y0 = lt(G_lower), x1 = lt(Vg_GP_mean), y1 = lt(G_upper))
  points(diag(G)~diag(Vg_GP_mean), col = "tomato3", pch = 19)
  abline(0, 1, col = "blue")
  text(0.2, 0.12, "Identity", col = "blue")
  #abline(lm(lt(G)~lt(Vg_GP_mean)))
  abline(v = 0)
  abline(h = 0)
  text(0.03, 0.35, "Variances", col = "tomato3")
  text(0.025, -0.05, "Co-variances")
}; g_predict_GP()

line = 2.5
cex = 2
las =2
padj = -10.5
png("data/growth_cov_prediction_composite_GP.png", width = 900, height = 900)
layout(matrix(c(1,2,3, 4), 2, 2, byrow = TRUE))
corrplot.mixed(cov2cor(Va_GP_mean), upper = "ellipse", mar = c(0, 0, 0, 0))
mtext("A", side=2, line=line, cex=cex, las=las, padj = padj)
corrplot.mixed(cov2cor(Vd_GP_mean), upper = "ellipse", mar = c(0, 0, 0, 0))
mtext("B", side=2, line=line, cex=cex, las=las, padj = padj)
corrplot.mixed(cov2cor(Vg_GP_mean), upper = "ellipse", mar = c(0, 0, 0, 0))
mtext("C", side=2, line=line, cex=cex, las=las, padj = padj)
par(mar = c(4, 5, 1, 1))
g_predict_GP()
mtext("D", side=2, line=line, cex=cex, las=las, padj = padj)
dev.off()


### Regressions with selection and divergence

a_corrs = data.frame(marker = 1:32, 
                     betaCorr = abs(apply(a_effect_matrix[,growth_traits], 1, function(x) vectorCor(x * growth_phen_sd, beta))),
                     dzCorr = abs(apply(a_effect_matrix[,growth_traits], 1, function(x) vectorCor(x * growth_phen_sd, d_z))),
                     norm = apply(a_effect_matrix[,growth_traits], 1, function(x) Norm(x * growth_phen_sd)))
d_corrs = data.frame(marker = 1:32,
                     betaCorr = abs(apply(d_effect_matrix[,growth_traits], 1, function(x) vectorCor(x * growth_phen_sd, beta))),
                     dzCorr = abs(apply(d_effect_matrix[,growth_traits], 1, function(x) vectorCor(x * growth_phen_sd, d_z))),
                     norm = apply(d_effect_matrix[,growth_traits], 1, function(x) Norm(x * growth_phen_sd)))
write.csv(a_corrs, "./data/growth_additive_correlations_beta_dZ.csv")
write.csv(d_corrs, "./data/growth_dominance_correlations_beta_dZ.csv")
ggplot(crss, aes(class, value, fill = class)) + geom_violin()
additive_beta = 
  ggplot(a_corrs, aes(norm, betaCorr)) + 
  geom_point(size = 3) + 
  #geom_point(aes(color = Early_proportion), size = 3) + scale_color_viridis() + 
  geom_smooth(method = "lm", color = "black") + 
  labs(x = "Additive effect vector norm", y = expression(paste("Alligment with ", beta)))
dominance_beta = ggplot(d_corrs, aes(norm, betaCorr)) + geom_point(size = 3) + geom_smooth(method = "lm", color = "black") + 
  labs(x = "Dominance effect vector norm", y = expression(paste("Alligment with ", beta)))
additive_dz = ggplot(a_corrs, aes(norm, dzCorr)) + geom_point(size = 3) + geom_smooth(method = "lm", color = "black") + 
  labs(x = "Additive effect vector norm", y = expression("Alligment with divergence"))
dominance_dz = ggplot(d_corrs, aes(norm, dzCorr)) + geom_point(size = 3) + geom_smooth(method = "lm", color = "black") + 
  labs(x = "Dominance effect vector norm", y = expression("Alligment with divergence"))
regressions = plot_grid(additive_beta, dominance_beta, additive_dz, dominance_dz, labels = LETTERS[1:4])
save_plot("data/growth_effect_aligment_regressions.png", regressions, base_height = 4, base_aspect_ratio = 2, ncol = 2, nrow = 2)

lm(betaCorr~norm, data = a_corrs) %>% summary
lm(dzCorr~norm, data = a_corrs) %>% summary
lm(betaCorr~norm, data = d_corrs) %>% summary
lm(dzCorr~norm, data = d_corrs) %>% summary

#### Prediction

a_corrs = data.frame(marker = 1:353, 
                     betaCorr = abs(apply(a_effect_matrix_HC[,growth_traits], 1, 
                                          function(x) vectorCor(x * growth_phen_sd, beta_w))),
                     dzCorr = abs(apply(a_effect_matrix_HC[,growth_traits], 1, 
                                        function(x) vectorCor(x * growth_phen_sd, d_z))),
                     norm = apply(a_effect_matrix_HC[,growth_traits], 1, 
                                  function(x) Norm(x * growth_phen_sd)))
d_corrs = data.frame(marker = 1:353,
                     betaCorr = abs(apply(d_effect_matrix_HC[,growth_traits], 1, 
                                          function(x) vectorCor(x * growth_phen_sd, beta_w))),
                     dzCorr = abs(apply(d_effect_matrix_HC[,growth_traits], 1, 
                                        function(x) vectorCor(x * growth_phen_sd, d_z))),
                     norm = apply(d_effect_matrix_HC[,growth_traits], 1, 
                                  function(x) Norm(x * growth_phen_sd)))

additive_beta = 
  ggplot(a_corrs, aes(norm, betaCorr)) + 
  geom_point(size = 3) + 
  geom_smooth(method = "lm", color = "black") + 
  labs(x = "Additive effect vector norm", y = expression(paste("Alligment with ", beta)))
dominance_beta = ggplot(d_corrs, aes(norm, betaCorr)) + geom_point(size = 3) + geom_smooth(method = "lm", color = "black") + 
  labs(x = "Dominance effect vector norm", y = expression(paste("Alligment with ", beta)))
additive_dz = ggplot(a_corrs, aes(norm, dzCorr)) + geom_point(size = 3) + geom_smooth(method = "lm", color = "black") + 
  labs(x = "Additive effect vector norm", y = expression("Alligment with divergence"))
dominance_dz = ggplot(d_corrs, aes(norm, dzCorr)) + geom_point(size = 3) + geom_smooth(method = "lm", color = "black") + 
  labs(x = "Dominance effect vector norm", y = expression("Alligment with divergence"))
regressions = plot_grid(additive_beta, dominance_beta, additive_dz, dominance_dz, labels = LETTERS[1:4])
save_plot("data/growth_effect_aligment_regressions_GP.png", regressions, base_height = 4, base_aspect_ratio = 2, ncol = 2, nrow = 2)

lm(betaCorr~norm, data = a_corrs) %>% summary
lm(dzCorr~norm, data = a_corrs) %>% summary
lm(betaCorr~norm, data = d_corrs) %>% summary
lm(dzCorr~norm, data = d_corrs) %>% summary

### Predictions

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

library(scales)

growth_observed_plot = ggplot() + scale_x_discrete(labels = paste("Week", 1:7)) + labs(y = "Weekly growth (g)", x = "Start week") + geom_line(size = 1, data = filter(growth_prediction, Type == "Observed"), aes(trait, value, group = Line, color = Line)) + theme_cowplot()
save_plot("data/growth_LG_SM_F3.png", growth_observed_plot, base_height = 7, base_aspect_ratio = 2)

par(mar = c(0, 0, 0, 0))
g_plot = ~corrplot.mixed(cov2cor(G),       upper = "ellipse", mar = c(0, 0, 0, 0))
growth_curves_cov = plot_grid(growth_observed_plot, g_plot, labels = LETTERS[1:2], ncol = 2)
save_plot("data/growth_LG_SM_F3_covF3.png", growth_curves_cov, base_height = 6, base_aspect_ratio = 1.5, ncol = 2)

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


growth_pred_plot_SUR = ggplot() + 
  scale_x_discrete(labels = paste("Week", 1:7)) + 
  labs(y = "Weekly growth (g)", x = "Start week") + 
  geom_line(data = post, color = "gray", size = 0.5, 
            linetype = 1, alpha = 0.1, aes(trait, value, group = interaction(Type, Line, iterations))) + 
  geom_line(size = 1, data = growth_prediction, 
            aes(trait, value, group = interaction(Type, Line), color = Line, linetype = Type)) +
  theme(legend.position = "none")  + 
  scale_y_continuous(limits = c(0,9.5)) + 
  annotate("text", label = "QTL\n mapping", x= 1.1, y = 8)

w_ad = rstan::extract(full_HCp, pars = "w_ad")[[1]]
effect_matrix = aaply(w_ad, c(2, 3), mean)
SM_e = rowSums(-1 * as.matrix(effect_matrix)) * sapply(growth_phen[,growth_traits], sd) + 
    sapply(growth_phen[,growth_traits], mean)
LG_e = rowSums(as.matrix(effect_matrix)) * sapply(growth_phen[,growth_traits], sd) + 
    sapply(growth_phen[,growth_traits], mean)
post_HC = adply(w_ad, 1, posterior_predict)
growth_prediction_HC = reshape2::melt(data.frame(trait = as.factor(growth_traits), 
                                              SM_QTL = SM_e,
                                              SM_Observed = SM,
                                              LG_QTL = LG_e,
                                              F3_Observed = F3,
                                              LG_Observed = LG)) %>% separate(variable, c("Line", "Type"))
(growth_pred_plot_HC_full = ggplot() + 
  scale_x_discrete(labels = paste("Week", 1:7)) + 
  labs(y = "Weekly growth (g)", x = "Start week") + 
  geom_line(data = post_HC, color = "gray", size = 0.5, 
            linetype = 1, alpha = 0.1, aes(trait, value, group = interaction(Type, Line, iterations))) + 
  geom_line(size = 1, data = growth_prediction_HC, 
            aes(trait, value, group = interaction(Line, Type), color = Line, linetype = Type)) + 
  theme(legend.position = c(0.7, 0.7), legend.title = element_blank()) + 
  scale_linetype_discrete(labels = c("Observed", "Predicted")) +
  scale_color_discrete(labels = c("F3", "LG/J", "SM/J")) + 
  scale_y_continuous(limits = c(0,9.5)) + 
  annotate("text", label = "Genome\n prediction", x= 1.1, y = 8))

predictions = plot_grid(growth_pred_plot_SUR, growth_pred_plot_HC_full, labels = LETTERS[1:2])
save_plot("data/growth_LG_SM_F3_predictions.png", predictions, base_height = 6, base_aspect_ratio = 1, ncol = 2)


## Pleiotropic partition

a_pleiotropic_partition = calcPleitropicPartition(a_effect_matrix)

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
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = c(0.45, 0.9)) + 
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

markerMatrix$id = 1:353
significantMarkerPos = inner_join(significantMarkerMatrix, markerMatrix, by = c("chrom", "marker"))

pleiotropic_partition = a_effect_matrix_HC
pleiotropic_partition[growth_traits] = sqrt(pleiotropic_partition[growth_traits]^2)
pleiotropic_partition$id = factor(pleiotropic_partition$id, levels = pleiotropic_partition$id)
a_pleiotropic_partition = pleiotropic_partition %>% 
  separate(id, c("chrom", "marker")) %>% 
  mutate(chrom = as.numeric(chrom), 
         marker = as.numeric(marker)) %>%
  arrange(chrom, marker) %>%
  mutate(id = factor(paste(chrom, marker, sep= "_"), levels = paste(chrom, marker, sep= "_")))

markerMatrix$id = 1:353
significantMarkerPos = inner_join(significantMarkerMatrix, markerMatrix, by = c("chrom", "marker"))
significantMarkerPos$size_ad = rowSums(a_pleiotropic_partition[significantMarkerPos$id, growth_traits])

a_pleiotropic_partition_plot =  
  ggplot() + 
  geom_bar(data = gather(a_pleiotropic_partition, key, value, growth_traits), 
           aes(id, value, group = key, color = key, fill = key), stat = "identity") +
  scale_y_continuous(limits = c(0, 0.35)) + background_grid(major = "x", minor = "none") + 
  scale_fill_viridis(discrete = TRUE, option = "D", 
                     guide = guide_legend(direction = "horizontal", 
                                          label.position = "top",
                                          nrow = 1, title = NULL)) + 
  scale_color_viridis(discrete = TRUE, option = "D", 
                      guide = guide_legend(direction = "horizontal",
                                           title = NULL)) + 
  geom_point(data = significantMarkerPos, aes(x = id, y = size_ad + 0.015), color = "tomato3", size = 2) +
  theme(legend.position = c(0.45, 0.9)) + 
  labs(x = "Chromossome Start", y = "Squared contribution to\n scaled pleitropic vector") +
  scale_x_discrete(breaks = a_pleiotropic_partition[markerMatrix[markerMatrix$marker==1,"id"],"id"],
                   labels = 1:19)
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

significantMarkerPos$size_dm = rowSums(d_pleiotropic_partition[significantMarkerPos$id, growth_traits])


d_pleiotropic_partition_plot =  
  ggplot() + 
  geom_bar(data = gather(d_pleiotropic_partition, key, value, growth_traits), aes(id, value, group = key, color = key, fill = key), stat = "identity") +
  scale_fill_viridis(discrete = TRUE, option = "D") + scale_color_viridis(discrete = TRUE, option = "D") + 
  scale_y_continuous(limits = c(0, 0.35)) + background_grid(major = "x", minor = "none") + 
  theme(legend.position = "none") +
  geom_point(data = significantMarkerPos, aes(x = id, y = size_dm + 0.015), color = "tomato3", size = 2) +
  labs(x = "Marker", y = "Squared contribution to\n scaled pleitropic vector") +
  scale_x_discrete(breaks = a_pleiotropic_partition[markerMatrix[markerMatrix$marker==1,"id"],"id"],
                   labels = 1:19)
pleiotropic_Effects_ad_dm = plot_grid(a_pleiotropic_partition_plot, d_pleiotropic_partition_plot, ncol = 1, labels = c("A", "B"))
pleiotropic_Effects_ad_dm = plot_grid(a_pleiotropic_partition_plot, d_pleiotropic_partition_plot, ad_biplot_HC, ncol = 1, labels = c("A", "B"))
save_plot("data/growth_pleiotropic_partition_ad_dm_GP.png", pleiotropic_Effects_ad_dm, base_height = 4.5, base_aspect_ratio = 2.5, nrow = 3) 
#save_plot("data/growth_pleiotropic_partition_dominance.png", pleiotropic_partition_plot, base_height = 7, base_aspect_ratio = 2) 


# Vizualising pleitropy per marker
load(file = "./Rdatas/significant_stan_fit.Rdata")
getStanEffects = function(markerMatrix, stan_model, trait_vector,
                          J = nrow(markerMatrix),
                          markers = 1:loci_per_chrom[current_chrom])
{
  K = length(trait_vector)
  HC_summary = summary(stan_model, pairs = c("w_ad", "w_dm"))$summary
  s = J * K * 2
  mask = grepl("w_", rownames(HC_summary))
  effects = data.frame(HC_summary[mask, c("mean", "2.5%", "97.5%")])
  colnames(effects) <- c("mean", "lower", "upper")
  effects$class = rep(c("additive", "dominance"), each = s/2)
  effects$chrom = rep(markerMatrix$chrom, K)
  effects$marker = rep(markerMatrix$marker, 2*K)
  effects$trait = rep(trait_vector, each = J)
  effects$id = factor(paste(effects$chrom, effects$marker, sep= "_"), 
                      levels = unique(paste(effects$chrom, effects$marker, sep= "_")))
  tbl_df(effects)
}
effectsStan = getStanEffects(significantMarkerMatrix, stan_model_SUR, growth_traits)
(effectsStan_plot = ggplot(effectsStan, aes(trait, mean, group = class, color = class)) + 
    geom_point(position = position_dodge(width = 0.2)) + 
    geom_pointrange(aes(ymin = lower, ymax = upper), position = position_dodge(width = 0.2)) + 
    geom_hline(yintercept = 0) + 
    scale_color_manual(values = c("black", "tomato3")) +
    facet_wrap(~id, scale = "free_y", ncol = 4) + 
  labs(x = "Week", y = "QTL effect") + scale_x_discrete(labels = 1:7))
save_plot("./data/growth_per_marker_additive_dominance_vectors_QTL.png", effectsStan_plot, base_height = 3, base_aspect_ratio = 1.3, 
          ncol = 4, nrow = 5)



## Comparison of QTL effects with GP effects for the same loci
a_effects_QTL_m = 
  gather(a_effect_matrix, trait, value, growth_traits) %>% 
  mutate(type = "QTL") %>%
  arrange(id, trait)
a_effects_GP_m = 
  gather(a_effect_matrix_HC, trait, value, growth_traits) %>% 
  mutate(type = "GP") %>%
  filter(id %in% a_effects_QTL_m$id) %>% 
  arrange(id, trait)
par(mfrow = c(1, 1), mar = c(4, 4, 4, 4))
plot(a_effects_QTL_m$value, a_effects_GP_m$value, pch = 19)
abline(0, 1)


## Comparison of direction between additive and dominance effects
i = 1

hist(laply(1:32, function(i) abs(vectorCor(as.numeric(a_effect_matrix[i,growth_traits]), 
                                       as.numeric(d_effect_matrix[i,growth_traits]))[1])), breaks = 30,
     main = "Additive-Dominance correlation distribution", xlab = "Correlations")

