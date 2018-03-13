setwd("/home/diogro/projects/mouse-qtls")
source("stanFunctions.R")
source('read_mouse_data.R')
source('utils.R')

stan_input = getStanInputFullGenome(growth_data, growth_traits)
names(stan_input)

full_HCp = readRDS("./Rdatas/growth_scaled_allmarkers_HCPlus")
str(full_HCp)

growth_data = inner_join(growth_phen_std, growth_markers, by = "ID")
growth_data[growth_traits] = scale(growth_data[growth_traits])

w0 = colMeans(rstan::extract(full_HCp, "w0")[[1]])
beta_family = colMeans(rstan::extract(full_HCp, "beta_family")[[1]])
w_ad = colMeans(rstan::extract(full_HCp, "w_ad")[[1]])
w_dm = colMeans(rstan::extract(full_HCp, "w_dm")[[1]])

w0m = matrix(w0, length(stan_input$family), 7, byrow = TRUE)

growth_sds = apply(growth_phen[,growth_traits], 2, sd)
G_HCp = cov(beta_family[stan_input$family,] + stan_input$ad %*% t(w_ad) + stan_input$dm %*% t(w_dm)) * outer(growth_sds, growth_sds)
dimnames(G_HCp) = NULL

old.par = par()
png("./data/growth_genomePrediction_fullSib_Gcorrelation.png", width = 1500, height = 800)
par(mfrow = c(1, 2), cex=2, oma = c(0, 0, 0, 0))
corrplot.mixed(cov2cor(G),       upper = "ellipse", mar = c(0, 0, 1, 0), main = "Family Full-Sib G")
corrplot.mixed(cov2cor(G_HCp), upper = "ellipse", mar = c(0, 0, 1, 0), main = "Genome prediction")
par(old.par)
dev.off()

MatrixCompare(G_HCp, G_stan)
