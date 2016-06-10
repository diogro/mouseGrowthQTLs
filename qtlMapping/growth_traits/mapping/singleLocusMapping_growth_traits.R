setwd("/home/diogro/projects/mouse-qtls")
source('read_mouse_data.R')
source('OAuth_lem_server.R')
1

Rdatas_folder = "~/gdrive/LGSM_project_Rdatas/"
#Rdatas_folder = "./data/Rdatas/"

install_load("MCMCglmm","doMC")
registerDoMC(40)

growth_data = inner_join(growth_phen_std,
                         Reduce(inner_join, growth_markers),
                         by = "ID")

value = paste("cbind(", paste(growth_traits, collapse = ', '), ")", sep = '')

fixed_effects = "trait - 1"

null_formula = paste(value, fixed_effects, sep = ' ~ ')

runNullMCMCModel <- function(null_formula, pl = TRUE, ...) {
    prior = list(R = list(V = diag(num_growth_traits), n = 0.002),
                 G = list(G1 = list(V = diag(num_growth_traits) * 0.02, n = 0.001)))
    growth_MCMC_null_model = MCMCglmm(as.formula(null_formula),
                                random = ~us(trait):FAMILY,
                                data = as.data.frame(growth_data),
                                rcov = ~us(trait):units,
                                family = rep("gaussian", num_growth_traits),
                                prior = prior,
                                pl = pl,
                                verbose = TRUE, ...)
    return(growth_MCMC_null_model)
}

#growth_MCMC_null_model = runNullMCMCModel(null_formula, nitt=150000, thin=100, burnin=50000)
#save(growth_MCMC_null_model, file = paste0(Rdatas_folder, "growth_MCMC_null_model.Rdata"))
load(paste0(Rdatas_folder, "growth_MCMC_null_model.Rdata"))
summary(growth_MCMC_null_model)

Gs_mcmc = array(growth_MCMC_null_model$VCV[,1:(num_growth_traits*num_growth_traits)], dim = c(1000, num_growth_traits, num_growth_traits))
G_mcmc = apply(array(growth_MCMC_null_model$VCV[,1:(num_growth_traits*num_growth_traits)], dim = c(1000, num_growth_traits, num_growth_traits)), 2:3, median)

R_mcmc = apply(array(growth_MCMC_null_model$VCV[,-c(1:(num_growth_traits*num_growth_traits))],
                     dim = c(1000, num_growth_traits, num_growth_traits)), 2:3, median)

makeMarkerList = function(pos) paste(paste('trait:chrom', pos[1],"_", c('A', 'D'), pos[2], sep = ''), collapse = ' + ')
markerMatrix = ldply(1:19, function(x) data.frame(chrom = x, locus = 1:loci_per_chrom[[x]]))
markerList = alply(markerMatrix, 1, makeMarkerList)

runSingleLocusMCMCModel <- function(marker_term, null_formula, start = NULL, ...){
    genotype.formula = paste(null_formula, marker_term, sep = ' + ')
    prior = list(R = list(V = diag(num_growth_traits), n = 0.002),
                 G = list(G1 = list(V = diag(num_growth_traits) * 0.02, n = num_growth_traits+1)))
    growth_MCMC_singleLocus = MCMCglmm(as.formula(genotype.formula),
                                random = ~us(trait):FAMILY,
                                data = as.data.frame(growth_data),
                                rcov = ~us(trait):units,
                                family = rep("gaussian", num_growth_traits),
                                start = start,
                                prior = prior,
                                verbose = FALSE,
                                ...)
    return(growth_MCMC_singleLocus)
}

start <- list(R = list(V = R_mcmc), G = list(G1 = G_mcmc), liab = matrix(growth_MCMC_null_model$Liab[1,], ncol = num_growth_traits))

all_loci_MCMC = llply(markerList, runSingleLocusMCMCModel, null_formula, start, nitt=153000, thin=10, burnin=3000, .parallel = TRUE)
save(all_loci_MCMC, file = paste0(Rdatas_folder, "growth_singleLocusMapping.Rdata"))
load(paste0(Rdatas_folder, "growth_singleLocusMapping.Rdata"))
x = all_loci_MCMC[[1]]
all_effects_mcmc = ldply(all_loci_MCMC,
      function(x){
           ad = summary(x)$solutions[(num_growth_traits+1):(num_growth_traits+num_growth_traits), ]
           colnames(ad) = paste("ad", colnames(ad), sep = "_")
           dm = summary(x)$solutions[(2*num_growth_traits+1):(2*num_growth_traits+num_growth_traits), ]
           colnames(dm) = paste("dm", colnames(dm), sep = "_")
           p_ad = colSums((x$Sol[,(num_growth_traits+1):(num_growth_traits+num_growth_traits)]) > 0)/nrow(x$Sol)
           p_dm = colSums((x$Sol[,(2*num_growth_traits+1):(2*num_growth_traits+num_growth_traits)]) > 0)/nrow(x$Sol)
           pos = na.omit(as.numeric(unlist(strsplit(unlist(as.character(x$Fixed$formula)[3]), "[^0-9]+"))))[2:3]
           data.frame(chrom = pos[1], marker = pos[2], trait = 1:num_growth_traits, ad, dm, p_ad, p_dm)
      }, .id = NULL, .parallel = TRUE) %>% tbl_df %>%
      dplyr::mutate(count = rep(seq(markerList), each = num_growth_traits)) %>%
      dplyr::select(count, everything()) %>% dplyr::select(-locus)
write_csv(all_effects_mcmc, "./data/growth traits/effectsSingleLocus.csv")
all_effects_mcmc = read_csv("./data/growth traits/effectsSingleLocus.csv")
dmSend("Finished growth single locus mapping", "diogro")

singleLocus_DIC_df <- 
  ldply(all_loci_MCMC, `[[`, "DIC") %>% 
  mutate(singleLocus_DICDiff = growth_MCMC_null_model$DIC - V1) %>% 
  rename(singleLocus_DIC = V1) %>% tbl_df
write_csv(singleLocus_DIC_df, "./data/growth traits/singleLocusDIC.csv")