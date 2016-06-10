setwd("/home/diogro/projects/mouse-qtls")
source('read_mouse_data.R')
source('OAuth_lem_server.R')
1

install_load("MCMCglmm","doMC")
registerDoMC(40)

#Rdatas_folder = "~/gdrive/LGSM_project_Rdatas/"
Rdatas_folder = "./data/Rdatas/"

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

G_mcmc = apply(array(growth_MCMC_null_model$VCV[,1:(num_growth_traits*num_growth_traits)], dim = c(1000, num_growth_traits, num_growth_traits)), 2:3, median)
R_mcmc = apply(array(growth_MCMC_null_model$VCV[,-c(1:(num_growth_traits*num_growth_traits))], dim = c(1000, num_growth_traits, num_growth_traits)), 2:3, median)

simulated_DIC_array = array(NA, c(1000, 20))
sim_chrom_number = 1
for(sim_chrom_number in 1:20){
    print(sim_chrom_number)
    ptm <- proc.time()
    growth_data = inner_join(growth_phen_std, simulated_markers[[sim_chrom_number]], by = "ID")

    makeMarkerList = function(pos) paste(paste('trait:sim_chrom', pos[1],"_", c('A', 'D'), pos[2], sep = ''), collapse = ' + ')
    markerMatrix = ldply(sim_chrom_number, function(x) data.frame(chrom = x, locus = 1:simulated_loci_per_chrom[[x]]))
    markerList = alply(markerMatrix, 1, makeMarkerList)

    runSingleLocusMCMCModelDIC <- function(marker_term, null_formula, start = NULL, verbose = FALSE, ...){
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
                                         verbose = verbose,
                                         ...)
        return(growth_MCMC_singleLocus$DIC)
    }

    start <- list(R = list(V = R_mcmc), G = list(G1 = G_mcmc), liab = matrix(growth_MCMC_null_model$Liab[1,], ncol = num_growth_traits))

    simulated_DIC = laply(markerList, runSingleLocusMCMCModelDIC, null_formula, start, nitt=3300, thin=15, burnin=300, .parallel = TRUE)
    simulated_DIC_array[,sim_chrom_number] = simulated_DIC
    time = proc.time() - ptm
    dmSend(paste("Finished simulated chromossome", sim_chrom_number, "in", round(time[3]/60, 2), "minutes." ), "diogro")
}

growth_nullDIC = growth_MCMC_null_model$DIC - simulated_DIC_array
saveRDS(growth_nullDIC, file = "./data/growth traits/nullDIC.rds")
dmSend("Finished null analysis", "diogro")

