setwd("~/projects/mouse-qtls/")
source("~/projects/mouse-qtls/read_mouse_data.R")
source("~/projects/mouse-qtls/OAuth_lem_server.R")

install_load("MCMCglmm","doMC")
registerDoMC(3)


makeSimData = function(data, marker, percent){
  a = findA(data, marker, percent)
  data + marker * a
}
markerVar = function(a, data, marker){
    data = mean(data) + residuals(lm(data ~ marker))
    za = data + marker * a
    n = length(data)
    V_pa = sum((za - mean(za))^2)/(n - 1)
    gen_means = tapply(za, marker, mean)
    gen_n = table(marker) / (n - 3)
    var_marker = gen_n %*% (gen_means - mean(za))^2
    (var_marker / V_pa)[1]
}
findA = function(data, marker, percent){
  V_p = var(data)
  a = sqrt((2*percent*V_p)/(1 - percent))
  return(a)
}

# data = rnorm(1000, 1, 10)
# marker = sample(c(-1, 0, 1), 1000, replace = TRUE)
# a = findA(data, marker, 0.01)
# markerVar(a, data, marker)

# sim_chrom_number = 1
# locus = 1
# effect_size = 0.01

simulateDIC = function(locus, sim_chrom_number, effect_size){
  print
  area_data = inner_join(area_phen_std, simulated_markers[[sim_chrom_number]], by = "ID")

  marker_column_A = data.frame(select(area_data, matches(paste0("_A", locus, "$"))))[,1]
  marker_column_D = data.frame(select(area_data, matches(paste0("_D", locus, "$"))))[,1]
  sim_data = tbl_df(as.data.frame(llply(select(area_data, area1:area7),
                                        makeSimData, marker_column_A, effect_size))) %>%
    mutate(A = marker_column_A, D = marker_column_D, FAMILY = area_phen_std$FAMILY)

  print(paste(sim_chrom_number, locus))
  value = paste("cbind(", paste(area_traits, collapse = ', '), ")", sep = '')

  fixed_effects = "trait - 1"

  null_formula = paste(value, fixed_effects, sep = ' ~ ')

  prior = list(R = list(V = diag(num_area_traits), n = 0.002),
               G = list(G1 = list(V = diag(num_area_traits) * 0.02, n = 0.001)))
  null_model = MCMCglmm(as.formula(null_formula),
                        random = ~us(trait):FAMILY,
                        data = as.data.frame(sim_data),
                        rcov = ~us(trait):units,
                        family = rep("gaussian", num_area_traits),
                        prior = prior,
                        nitt=8000, thin=20, burnin=3000,
                        pl = TRUE,
                        verbose = FALSE)
  n_mc = dim(null_model$Sol)[1]
  G_mcmc = apply(array(null_model$VCV[,1:(num_area_traits*num_area_traits)], dim = c(n_mc, num_area_traits, num_area_traits)), 2:3, median)
  R_mcmc = apply(array(null_model$VCV[,-c(1:(num_area_traits*num_area_traits))], dim = c(n_mc, num_area_traits, num_area_traits)), 2:3, median)
  start <- list(R = list(V = R_mcmc), G = list(G1 = G_mcmc), liab = matrix(null_model$Liab[1,], ncol = num_area_traits))

  genotype.formula = paste(null_formula, "trait:A + trait:D", sep = ' + ')
  prior = list(R = list(V = diag(num_area_traits), n = 0.002),
               G = list(G1 = list(V = diag(num_area_traits) * 0.02, n = num_area_traits+1)))
  marker_model = MCMCglmm(as.formula(genotype.formula),
                          random = ~us(trait):FAMILY,
                          data = as.data.frame(sim_data),
                          rcov = ~us(trait):units,
                          family = rep("gaussian", num_area_traits),
                          nitt=1300, thin=5, burnin=300,
                          start = start,
                          prior = prior,
                          verbose = FALSE)
  null_model$DIC - marker_model$DIC
}

effect_size = seq(0.005, 0.04, length.out = 20)
DIC_list = vector("list", 20)
names(DIC_list) = effect_size
for(i in 1:20){
    ptm <- proc.time()
    DIC_list[[i]] = ldply(1:3, simulateDIC, i, effect_size[i], .parallel = TRUE)
    time = proc.time() - ptm
    dmSend(paste("Finished simulated chromossome", i, "in", round(time[3]/60, 2), "minutes." ), "diogro")
}
DIC_power = ldply(DIC_list)
write_csv(DIC_power, "./data/area traits/power_analysis.csv")
dmSend("Finished power analysis", "diogro")
