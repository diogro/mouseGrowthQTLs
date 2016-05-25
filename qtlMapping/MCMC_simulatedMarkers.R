setwd("/home/diogro/projects/mouse-qtls")
source('read_mouse_data.R')

install_load("MCMCglmm","doMC")
registerDoMC(4)

sim_chrom_number = 2

area_data = inner_join(area_phen_std, simulated_markers[[sim_chrom_number]], by = "ID")

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

#area_MCMC_null_model = runNullMCMCModel(null_formula, nitt=150000, thin=100, burnin=50000)
#save(area_MCMC_null_model, file = "./data/Rdatas/area_MCMC_null_model.Rdata")
load("./data/Rdatas/area_MCMC_null_model.Rdata")
summary(area_MCMC_null_model)

Gs_mcmc = array(area_MCMC_null_model$VCV[,1:(num_area_traits*num_area_traits)], dim = c(1000, num_area_traits, num_area_traits))
G_mcmc = apply(array(area_MCMC_null_model$VCV[,1:(num_area_traits*num_area_traits)], dim = c(1000, num_area_traits, num_area_traits)), 2:3, median)

R_mcmc = apply(array(area_MCMC_null_model$VCV[,-c(1:(num_area_traits*num_area_traits))],
                     dim = c(1000, num_area_traits, num_area_traits)), 2:3, median)

makeMarkerList = function(pos) paste(paste('trait:sim_chrom', pos[1],"_", c('A', 'D'), pos[2], sep = ''), collapse = ' + ')
markerMatrix = ldply(2, function(x) data.frame(chrom = x, locus = 1:simulated_loci_per_chrom[[x]]))
markerList = alply(markerMatrix, 1, makeMarkerList)
markerList[[1]]

runSingleLocusMCMCModel <- function(marker_term, null_formula, start = NULL, verbose = FALSE, ...){
    genotype.formula = paste(null_formula, marker_term, sep = ' + ')
    prior = list(R = list(V = diag(num_area_traits), n = 0.002),
                 G = list(G1 = list(V = diag(num_area_traits) * 0.02, n = num_area_traits+1)))
    area_MCMC_singleLocus = MCMCglmm(as.formula(genotype.formula),
                                random = ~us(trait):FAMILY,
                                data = as.data.frame(area_data),
                                rcov = ~us(trait):units,
                                family = rep("gaussian", num_area_traits),
                                start = start,
                                prior = prior,
                                verbose = verbose,
                                ...)
    return(area_MCMC_singleLocus)
}

start <- list(R = list(V = R_mcmc), G = list(G1 = G_mcmc), liab = matrix(area_MCMC_null_model$Liab[1,], ncol = num_area_traits))

all_loci_MCMC = llply(markerList, runSingleLocusMCMCModel, null_formula, start, nitt=3300, thin=15, burnin=300, .parallel = TRUE)
model_file = paste0("./data/Rdatas/simulated_mcmc_data", sim_chrom_number, ".Rdata")
save(all_loci_MCMC, file = model_file)
load(model_file)
x = all_loci_MCMC[[1]]
all_effects_mcmc = ldply(all_loci_MCMC,
      function(x){
           ad = summary(x)$solutions[8:14, ]
           colnames(ad) = paste("ad", colnames(ad), sep = "_")
           dm = summary(x)$solutions[15:21, ]
           colnames(dm) = paste("dm", colnames(dm), sep = "_")
           p_ad = colSums((x$Sol[,8:14]) > 0)/nrow(x$Sol)
           p_dm = colSums((x$Sol[,15:21]) > 0)/nrow(x$Sol)
           pos = na.omit(as.numeric(unlist(strsplit(unlist(as.character(x$Fixed$formula)[3]), "[^0-9]+"))))[2:3]
           data.frame(chrom = pos[1], marker = pos[2], trait = 1:num_area_traits, ad, dm, p_ad, p_dm)
      }, .id = NULL, .parallel = TRUE) %>% tbl_df %>%
      dplyr::mutate(count = rep(seq(markerList), each = 7)) %>%
      dplyr::select(count, everything()) %>% dplyr::select(-locus)
effect_file = paste0("./data/area traits/simulatedEffects_data", sim_chrom_number, "_mcmc.csv")
write_csv(all_effects_mcmc, effect_file)
all_effects_mcmc = read_csv(effect_file)
