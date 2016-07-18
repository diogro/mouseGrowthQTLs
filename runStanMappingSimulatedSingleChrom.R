setwd("/home/diogro/projects/mouse-qtls")
source('stanFunctions.R')
source('read_mouse_data.R')
source('utils.R')

Rdatas_folder = "~/gdrive/LGSM_project_Rdatas/"
#Rdatas_folder = "./data/Rdatas/"

options(mc.cores = 3)

weight_data = inner_join(weight_phen_std, simulated_genomes[[8]], by = "ID")

#weight_data = inner_join(weight_phen_std, markers, by = "ID")

weight_data[weight_traits] = scale(weight_data[weight_traits])
#weight_data[weight_traits] = scale(weight_data[weight_traits], scale = rep(sqrt(1/20), num_weight_traits))

#trait_vector = weight_traits
#n_effects = 10
#current_chrom = 6
#current_data = weight_data
#eff_sd = 0.1
simulateData = function(current_chrom, current_data, trait_vector,
                        n_effects = length(trait_vector), eff_sd = 0.1){
    true_effects = NULL
    for(i in 1:n_effects){
        current_trait = sample(trait_vector, 1)
        current_loci = sample(loci_per_chrom[current_chrom], 1)
        type = sample(c("A", "D"), 1)
        locus = paste0("chrom", current_chrom, "_", type, current_loci)
        effect = rnorm(1, 0, eff_sd)
        current_data[current_trait] = current_data[current_trait] + effect * current_data[locus]
        true_effects = rbind(true_effects, data.frame(chrom = current_chrom, marker = current_loci,
                   type = ifelse(type == "A", "additive", "dominance"),
                   true_effects = effect, trait = current_trait))
    }
    return(list(current_data, true_effects))
}
x = simulateData(current_chrom, weight_data, weight_traits, n_effects = 10)
weight_data = x[[1]]
true_effects = x[[2]]

weight_data[weight_traits] = scale(weight_data[weight_traits])
#weight_data[weight_traits] = scale(weight_data[weight_traits], scale = rep(sqrt(1/40), num_weight_traits))
apply(as.matrix(weight_data[weight_traits]), 2, var)

current_chrom = 6
sim_model = runStanModel(current_chrom, weight_data, weight_traits,
                         chain = 3, iter = 200, model_file = './SUR_horseShoePlus.stan')

plotEffectEstimate(current_chrom, sim_model[[1]], "weight")
plotWeights(current_chrom, sim_model[[2]], "weight")

names(sim_model[[3]])
colMeans(sim_model[[3]]$shrink_ad)
G = colMeans(sim_model[[3]]$G)
R = colMeans(sim_model[[3]]$R)
G + R
cov(weight_data[weight_traits])

plotEffectEstimate(current_chrom, sim_model[[1]], "sim_scaled_weight_0.1s_random_0.1", true_effects)
plotWeights(current_chrom, sim_model[[2]], "sim_weight_0.1s_random_0.1")
