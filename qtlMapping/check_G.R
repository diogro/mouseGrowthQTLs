install_load("evolqg")

G_sas = matrix(NA, 7, 7)
G_sas[1,1] = 0.004738
G_sas[2,1] = 0.00843
G_sas[2,2] = 0.1159
G_sas[3,1] = 0.005293
G_sas[3,2] = 0.04881
G_sas[3,3] = 0.06629
G_sas[4,1] = 0.01916
G_sas[4,2] = 0.08291
G_sas[4,3] = 0.105
G_sas[4,4] = 0.2773
G_sas[5,1] = 0.01407
G_sas[5,2] = 0.07885
G_sas[5,3] = 0.07617
G_sas[5,4] = 0.173
G_sas[5,5] = 0.1424
G_sas[6,1] = 0.002259
G_sas[6,2] = 0.009621
G_sas[6,3] = 0.008647
G_sas[6,4] = 0.0227
G_sas[6,5] = 0.0165
G_sas[6,6] = 0.004291
G_sas[7,1] = 0.01108
G_sas[7,2] = 0.04479
G_sas[7,3] = 0.03933
G_sas[7,4] = 0.1017
G_sas[7,5] = 0.08045
G_sas[7,6] = 0.01106
G_sas[7,7] = 0.05909

G_sas[upper.tri(G_sas)] <- t(G_sas)[upper.tri(G_sas)]
p_cov = ggplot(data.frame(sas = G_sas[lower.tri(G_sas, diag = TRUE)], ML = G_ML[lower.tri(G_ML, diag = TRUE)]),
       aes(sas, ML)) + geom_point() + geom_smooth(method = "lm") + geom_abline() + ggtitle("Covariances")
p_cor = ggplot(data.frame(sas = cov2cor(G_sas)[lower.tri(G_sas, diag = TRUE)], ML = cov2cor(G_ML)[lower.tri(G_ML, diag = TRUE)]),
       aes(sas, ML)) + geom_point() + geom_smooth(method = "lm") + geom_abline() + ggtitle("Correlations")

(G_comparions_plot = plot_grid(p_cov, p_cor))
print(G_comparions_plot)
#save_plot("./data/area traits/G_comparison.png", G_comparions_plot, ncol = 3, nrow = 2, base_height = 5)
