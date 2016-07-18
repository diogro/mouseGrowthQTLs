setwd("/home/diogro/projects/mouse-qtls")
eource('read_mouse_data.R')
source('OAuth_lem_server.R')
1

#Rdatas_folder = "~/gdrive/LGSM_project_Rdatas/"
Rdatas_folder = "./data/Rdatas/"

install_load("MCMCglmm","doMC")
registerDoMC(80)

area_data = inner_join(area_phen_std, area_markers, by = "ID")

value = paste("cbind(", paste(area_traits, collapse = ', '), ")", sep = '')

fixed_effects = "trait - 1"

null_formula = paste(value, fixed_effects, sep = ' ~ ')

runNullMCMCModel <- function(null_formula, pl = TRUE, ...) {
  prior = list(R = list(V = diag(num_area_traits), n = 0.002),
               G = list(G1 = list(V = diag(num_area_traits) * 0.02, n = 0.001)))
  area_MCMC_null_model = MCMCglmm(as.formula(null_formula),
                                      random = ~us(trait):FAMILY,
                                      data = as.data.frame(area_data),
                                      rcov = ~us(trait):units,
                                      family = rep("gaussian", num_area_traits),
                                      prior = prior,
                                      pl = pl,
                                      verbose = TRUE, ...)
  return(area_MCMC_null_model)
}

area_MCMC_null_model = runNullMCMCModel(null_formula, nitt=1003000, thin=1000, burnin=3000)
save(area_MCMC_null_model, file = paste0(Rdatas_folder, "area_MCMC_null_model.Rdata"))
load(paste0(Rdatas_folder, "area_MCMC_null_model.Rdata"))
summary(area_MCMC_null_model)

G_mcmc = apply(array(area_MCMC_null_model$VCV[,1:(num_area_traits*num_area_traits)], dim = c(1000, num_area_traits, num_area_traits)), 2:3, median)
R_mcmc = apply(array(area_MCMC_null_model$VCV[,-c(1:(num_area_traits*num_area_traits))], dim = c(1000, num_area_traits, num_area_traits)), 2:3, median)

source("./qtlMapping/area_traits/mapping/singleLocusMapping_area_traits.R")

flank_dist = 5

source("./qtlMapping/area_traits/mapping/intervalMapping_area_traits.R")

model_file = paste0(Rdatas_folder, "area_intervalMapping_", flank_dist, "cM_mcmc.Rdata")
load(model_file)

intervalMapping_DIC_df <- 
  ldply(intervalMapping_MCMC, function(x) c(flankDIC_5cM = x[[1]]$DIC, 
                                            focalDIC_5cM = x[[2]]$DIC)) %>% 
  mutate(DICDiff_5cM = flankDIC_5cM - focalDIC_5cM) %>% tbl_df
rm(intervalMapping_MCMC); gc()

flank_dist = 10

source("./qtlMapping/area_traits/mapping/intervalMapping_area_traits.R")

model_file = paste0(Rdatas_folder, "area_intervalMapping_", flank_dist, "cM_mcmc.Rdata")
load(model_file)

intervalMapping_DIC_df <- 
  ldply(intervalMapping_MCMC, function(x) c(flankDIC_10cM = x[[1]]$DIC, 
                                            focalDIC_10cM = x[[2]]$DIC)) %>% 
  mutate(DICDiff_10cM = flankDIC_10cM - focalDIC_10cM) %>% tbl_df %>% 
  inner_join(intervalMapping_DIC_df, by = c("chrom", "marker")) 

rm(intervalMapping_MCMC); gc()

flank_dist = 15

source("./qtlMapping/area_traits/mapping/intervalMapping_area_traits.R")

model_file = paste0(Rdatas_folder, "area_intervalMapping_", flank_dist, "cM_mcmc.Rdata")
load(model_file)

intervalMapping_DIC_df <- 
  ldply(intervalMapping_MCMC, function(x) c(flankDIC_15cM = x[[1]]$DIC, 
                                            focalDIC_15cM = x[[2]]$DIC)) %>% 
  mutate(DICDiff_15cM = flankDIC_15cM - focalDIC_15cM) %>% tbl_df %>% 
  inner_join(intervalMapping_DIC_df, by = c("chrom", "marker")) 

rm(intervalMapping_MCMC); gc()

flank_dist = 20

source("./qtlMapping/area_traits/mapping/intervalMapping_area_traits.R")

model_file = paste0(Rdatas_folder, "area_intervalMapping_", flank_dist, "cM_mcmc.Rdata")
load(model_file)

intervalMapping_DIC_df <- 
  ldply(intervalMapping_MCMC, function(x) c(flankDIC_20cM = x[[1]]$DIC, 
                                            focalDIC_20cM = x[[2]]$DIC)) %>% 
  mutate(DICDiff_20cM = flankDIC_20cM - focalDIC_20cM) %>% tbl_df %>% 
  inner_join(intervalMapping_DIC_df, by = c("chrom", "marker")) 

rm(intervalMapping_MCMC); gc()

write_csv(intervalMapping_DIC_df, "./data/area traits/intervalMappingDIC.csv")
