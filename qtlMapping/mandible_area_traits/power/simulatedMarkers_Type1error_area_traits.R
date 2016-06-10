setwd("/home/diogro/projects/mouse-qtls")
source('read_mouse_data.R')

install_load("MCMCglmm","doMC")
registerDoMC(40)

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
#save(area_MCMC_null_model, file = paste0(Rdatas_folder, "area_MCMC_null_model.Rdata"))
load(paste0(Rdatas_folder, "area_MCMC_null_model.Rdata"))
summary(area_MCMC_null_model)

G_mcmc = apply(array(area_MCMC_null_model$VCV[,1:(num_area_traits*num_area_traits)], dim = c(1000, num_area_traits, num_area_traits)), 2:3, median)
R_mcmc = apply(array(area_MCMC_null_model$VCV[,-c(1:(num_area_traits*num_area_traits))], dim = c(1000, num_area_traits, num_area_traits)), 2:3, median)

simulated_DIC_array = array(NA, c(1000, 20))
sim_chrom_number = 1
for(sim_chrom_number in 3:20){
    print(sim_chrom_number)
    area_data = inner_join(area_phen_std, simulated_markers[[sim_chrom_number]], by = "ID")

    makeMarkerList = function(pos) paste(paste('trait:sim_chrom', pos[1],"_", c('A', 'D'), pos[2], sep = ''), collapse = ' + ')
    markerMatrix = ldply(sim_chrom_number, function(x) data.frame(chrom = x, locus = 1:simulated_loci_per_chrom[[x]]))
    markerList = alply(markerMatrix, 1, makeMarkerList)

    runSingleLocusMCMCModelDIC <- function(marker_term, null_formula, start = NULL, verbose = FALSE, ...){
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
        return(area_MCMC_singleLocus$DIC)
    }

    start <- list(R = list(V = R_mcmc), G = list(G1 = G_mcmc), liab = matrix(area_MCMC_null_model$Liab[1,], ncol = num_area_traits))

    simulated_DIC = laply(markerList, runSingleLocusMCMCModelDIC, null_formula, start, nitt=3300, thin=15, burnin=300, .parallel = TRUE)
    simulated_DIC_array[,sim_chrom_number] = simulated_DIC
}

area_nullDIC = area_MCMC_null_model$DIC - simulated_DIC_array
saveRDS(area_nullDIC, file = "./data/area traits/nullDIC.rds")

