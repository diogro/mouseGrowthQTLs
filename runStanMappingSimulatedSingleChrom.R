setwd("/home/diogro/projects/mouse-qtls")
source('stanFunctions.R')
source('read_mouse_data.R')
source('utils.R')

Rdatas_folder = "~/gdrive/LGSM_project_Rdatas/"
#Rdatas_folder = "./data/Rdatas/"

rstan_options(auto_write = TRUE)
options(mc.cores = 4)

weight_data = inner_join(weight_phen_std, simulated_genomes[[8]], by = "ID")

#weight_data = inner_join(weight_phen_std, markers, by = "ID")

weight_data[weight_traits] = scale(weight_data[weight_traits])

trait_vector = weight_traits
n_effects = 10
current_chrom = 6
current_data = weight_data
eff_sd = 0.1
simulateData = function(current_chrom, current_data, trait_vector,
                        n_effects = length(trait_vector), eff_sd = 0.1){
    true_effects = NULL
    for(i in 1:n_effects){
        current_trait = sample(trait_vector, 3)
        current_loci = sample(loci_per_chrom[current_chrom], 1)
        type = sample(c("A", "D"), 1)
        locus = paste0("chrom", current_chrom, "_", type, current_loci)
        effect = rnorm(1, 0, eff_sd)
        for(trait in current_trait){
            current_data[trait] = current_data[trait] + effect * current_data[locus]
            true_effects = rbind(true_effects, data.frame(chrom = current_chrom, marker = current_loci,
                                                          type = ifelse(type == "A", "additive", "dominance"),
                                                          true_effects = effect, trait = trait))
        }
    }
    return(list(current_data, true_effects))
}
x = simulateData(current_chrom <- 6, weight_data, weight_traits, n_effects = 10)
sim_weight_data = x[[1]]
true_effects = x[[2]]

#weight_data[weight_traits] = scale(weight_data[weight_traits])
apply(as.matrix(sim_weight_data[weight_traits]), 2, var)

sim_model = runStanModel(current_chrom, sim_weight_data, weight_traits,
                         chain = 3, iter = 200, model_file = './SUR_HAL.stan')

colMeans(sim_model[[3]]$w_ad)

png("./data/figures/lp.png")
plot(sim_model[[3]]$'lp__')
dev.off()

#plotEffectEstimate(current_chrom, sim_model[[1]], "weight")
#plotShrinkage(current_chrom, sim_model[[2]], "weight")

plotEffectEstimate(current_chrom, sim_model[[1]], "sim_scaled_weight_HAL_random_0.1", true_effects, "free")
plotShrinkage(current_chrom, sim_model[[2]], "sim_scaled_weight_HAL_random_0.1")
