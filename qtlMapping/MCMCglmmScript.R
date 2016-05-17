source('read_mouse_data.R')

install_load("MCMCglmm")

area_data = inner_join(area_phen_std,
                       Reduce(inner_join, area_markers),
                       by = "ID")

value = paste("cbind(", paste(area_traits, collapse = ', '), ")", sep = '')

fixed_effects = "trait - 1"

null_formula = paste(value, fixed_effects, sep = ' ~ ')

runNullMCMCModel <- function(null_formula, pl = TRUE) {
    prior = list(R = list(V = diag(num_area_traits), n = 0.002),
                 G = list(G1 = list(V = diag(num_area_traits) * 0.02, n = 0.001)))
    area_MCMC_null_model = MCMCglmm(as.formula(null_formula),
                                random = ~us(trait):FAMILY,
                                data = as.data.frame(area_data),
                                rcov = ~us(trait):units,
                                family = rep("gaussian", num_area_traits),
                                prior = prior,
                                pl = pl,
                                verbose = TRUE)
    return(area_MCMC_null_model)
}

area_MCMC_null_model = runNullMCMCModel(null_formula)
summary(area_MCMC_null_model)

G_mcmc = apply(array(area_MCMC_null_model$VCV[,1:(num_area_traits*num_area_traits)], 
                     dim = c(1000, num_area_traits, num_area_traits)), 2:3, median)

R_mcmc = apply(array(area_MCMC_null_model$VCV[,-c(1:(num_area_traits*num_area_traits))], 
                     dim = c(1000, num_area_traits, num_area_traits)), 2:3, median)

makeMarkerList = function(pos) paste(paste('trait:chrom', pos[1],"_", c('A', 'D'), pos[2], sep = ''), collapse = ' + ')
markerMatrix = ldply(1:19, function(x) data.frame(chrom = x, locus = 1:loci_per_chrom[[x]]))
markerList = alply(markerMatrix, 1, makeMarkerList)

marker_term = markerList[[40]]
runSingleLocusMCMCModel <- function(marker_term, null_formula, start = NULL){
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
                                verbose = FALSE)
    return(area_MCMC_singleLocus)
}

start <- list(R = list(V = R_mcmc), G = list(G1 = G_mcmc), liab = matrix(area_MCMC_null_model$Liab[1,], ncol = num_area_traits))

all_loci_MCMC = llply(markerList, runSingleLocusMCMCModel, null_formula, start, .parallel = TRUE)
