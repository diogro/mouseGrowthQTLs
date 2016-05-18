setwd("/home/diogro/projects/mouse-qtls")
source('read_mouse_data.R')

install_load("MCMCglmm","doMC")
registerDoMC(100)

area_data = inner_join(area_phen_std,
                       Reduce(inner_join, area_markers),
                       by = "ID")

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

area_MCMC_null_model = runNullMCMCModel(null_formula, nitt=150000, thin=100, burnin=50000)
summary(area_MCMC_null_model)

G_mcmc = apply(array(area_MCMC_null_model$VCV[,1:(num_area_traits*num_area_traits)],
                     dim = c(1000, num_area_traits, num_area_traits)), 2:3, median)

R_mcmc = apply(array(area_MCMC_null_model$VCV[,-c(1:(num_area_traits*num_area_traits))],
                     dim = c(1000, num_area_traits, num_area_traits)), 2:3, median)

makeMarkerList = function(pos) paste(paste('trait:chrom', pos[1],"_", c('A', 'D'), pos[2], sep = ''), collapse = ' + ')
markerMatrix = ldply(1:19, function(x) data.frame(chrom = x, locus = 1:loci_per_chrom[[x]]))
markerList = alply(markerMatrix, 1, makeMarkerList)

runSingleLocusMCMCModel <- function(marker_term, null_formula, start = NULL, ...){
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
                                verbose = FALSE,
                                ...)
    return(area_MCMC_singleLocus)
}

start <- list(R = list(V = R_mcmc), G = list(G1 = G_mcmc), liab = matrix(area_MCMC_null_model$Liab[1,], ncol = num_area_traits))

all_loci_MCMC = llply(markerList, runSingleLocusMCMCModel, null_formula, start, nitt=150000, thin=100, burnin=50000, .parallel = TRUE)
save(all_loci_MCMC, file = "./data/Rdatas/all_loci_mcmc.Rdata")
all_effects_mcmc = ldply(all_loci_MCMC,
      function(x){
           ad = summary(x)$solutions[8:14, ]
           colnames(ad) = paste("ad", colnames(ad), sep = "_")
           dm = summary(x)$solutions[15:21, ]
           colnames(dm) = paste("dm", colnames(dm), sep = "_")
           pos = na.omit(as.numeric(unlist(strsplit(unlist(as.character(x$Fixed$formula)[3]), "[^0-9]+"))))[2:3]
           data.frame(chrom = pos[1], marker = pos[2], trait = 1:num_area_traits, ad, dm)
      }, .id = NULL, .parallel = TRUE) %>% tbl_df %>%
      mutate(count = rep(seq(markerList), each = 7)) %>% select(count, everything()) %>% select(-locus)
write_csv(all_effects_mcmc, "./data/area traits/effects_mcmc.csv")

marker_term = markerList[[339]]
area_MCMC_singleLocus = runSingleLocusMCMCModel(marker_term, null_formula, start, nitt=20000, thin=100, burnin=10000)
summary(area_MCMC_singleLocus)
