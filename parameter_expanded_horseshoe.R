setwd("/home/diogro/projects/mouse-qtls")
source('read_mouse_data.R')

Rdatas_folder = "~/gdrive/LGSM_project_Rdatas/"
#Rdatas_folder = "./data/Rdatas/"


install_load("MCMCglmm","doMC")
registerDoMC(4)

area_data = inner_join(area_phen_std,
                       simulated_genomes[[2]],
                       by = "ID")

area_data = area_data %>% mutate_each(funs(scale), matches('area'))

true_effects = c(findA(area_data$area1, area_data$chrom4_A10, 0.05),
                 findA(area_data$area2, area_data$chrom4_A10, 0.05),
                 findA(area_data$area3, area_data$chrom4_A13, 0.05),
                 findA(area_data$area4, area_data$chrom4_A14, 0.05),
                 findA(area_data$area5, area_data$chrom4_A14, 0.05),
                 findA(area_data$area6, area_data$chrom4_A13, 0.05),
                 findA(area_data$area7, area_data$chrom4_A14, 0.05))
true_effects = data.frame(true_effects, trait = area_traits, type = "additive",
                          marker = c(rep(10, 2), 13, rep(14, 2), 13, 14))

area_data$area1 = makeSimData(area_data$area1, area_data$chrom4_A10, 0.05)
area_data$area2 = makeSimData(area_data$area2, area_data$chrom4_A10, 0.05)
area_data$area3 = makeSimData(area_data$area3, area_data$chrom4_A13, 0.05)
area_data$area4 = makeSimData(area_data$area4, area_data$chrom4_A14, 0.05)
area_data$area5 = makeSimData(area_data$area5, area_data$chrom4_A14, 0.05)
area_data$area6 = makeSimData(area_data$area6, area_data$chrom4_A13, 0.05)
area_data$area7 = makeSimData(area_data$area7, area_data$chrom4_A14, 0.05)


value = paste("cbind(", paste(area_traits, collapse = ', '), ")", sep = '')

fixed_effects = "trait - 1"

null_formula = paste(value, fixed_effects, sep = ' ~ ')

current_chrom = 4
random_effects = paste0("~us(trait):FAMILY + idh(",
                        paste(paste0(c(paste0("chrom", current_chrom, "_A"),
                                       paste0("chrom", current_chrom, "_D")),
                                     rep(1:loci_per_chrom[current_chrom], each = 2)),
                              collapse = " + "),
                        "):trait")

prior = list(R = list(V = diag(num_area_traits), n = 0.002),
             G = list(G1 = list(V = diag(num_area_traits) * 0.02, n = 0.001),
                      G2 = list(V = diag(2*loci_per_chrom[current_chrom])*0.01, nu = 1,
                                alpha.mu = rep(0, 2*loci_per_chrom[current_chrom]),
                                alpha.V = diag(2*loci_per_chrom[current_chrom]))))
area_HC_model = MCMCglmm(as.formula(null_formula),
                         random = as.formula(random_effects),
                         data = as.data.frame(area_data),
                         rcov = ~us(trait):units,
                         family = rep("gaussian", num_area_traits),
                         prior = prior,
                         nitt = 2000, burnin = 1000, thin = 10,
                         pr = TRUE,
                         verbose = TRUE)

#write_rds(area_HC_model, paste0(Rdatas_folder, "HC_chrom4_area_s100.rds"))
#write_rds(area_HC_model, paste0(Rdatas_folder, "HC_chrom4_area_s10.rds"))
#write_rds(area_HC_model, paste0(Rdatas_folder, "HC_chrom4_area_s1.rds"))

hc_models = list(
s100 = read_rds(paste0(Rdatas_folder, "HC_chrom4_area_s100.rds")),
s10 = read_rds(paste0(Rdatas_folder, "HC_chrom4_area_s10.rds")),
s1 = read_rds(paste0(Rdatas_folder, "HC_chrom4_area_s1.rds"))
)
getHCEffects = function(area_HC_model){
  HC_summary = summary(area_HC_model, random=TRUE)
  n = dim(HC_summary$solutions)[1]
  s = n - 2*loci_per_chrom[current_chrom]*num_area_traits +1
  effects = data.frame(HC_summary$solutions[s:n,])
  colnames(effects) <- c("mean", "lower", "upper", "eff_samp", "pMCMC")
  effects$type = rep(c("additive", "dominance"), each = num_area_traits)
  effects$chrom = current_chrom
  effects$marker = rep(1:loci_per_chrom[current_chrom], each = 2*num_area_traits)
  effects$trait = area_traits
  tbl_df(effects)
}
hc_effects = ldply(hc_models, getHCEffects)
effects = getHCEffects(area_HC_model)
hc_plot = ggplot(effects, aes(marker, mean, group = trait)) +
  geom_point() + facet_grid(trait~type, scales = "free") +
  geom_hline(yintercept = 0) +
  geom_point(size = 0.3) +
  geom_point(data = true_effects, aes(y = true_effects), color = "red") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, size = 0.3)
save_plot("data/figures/sim_sparce_regression_s100_e0.03.png", hc_plot, base_height = 6, base_aspect_ratio = 1.8)
