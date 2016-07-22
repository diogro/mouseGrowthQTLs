setwd("/home/diogro/projects/mouse-qtls")
source('stanFunctions.R')
source('read_mouse_data.R')
source('utils.R')

Rdatas_folder = "~/gdrive/LGSM_project_Rdatas/"
#Rdatas_folder = "./data/Rdatas/"

rstan_options(auto_write = TRUE)
options(mc.cores = 3)

growth_data = inner_join(growth_phen_std, simulated_genomes[[8]], by = "ID")

growth_data = inner_join(growth_phen_std, markers, by = "ID")

growth_data[growth_traits] = scale(growth_data[growth_traits])
simulateData = function(current_chrom, current_data, trait_vector,
                        n_effects = length(trait_vector), eff_sd = 0.1){
    true_effects = NULL
    all_loci = sample(loci_per_chrom[current_chrom], n_effects)
    for(i in 1:n_effects){
        current_trait = sample(trait_vector, 3)
        type = sample(c("A", "D"), 1)
        current_loci = all_loci[i]
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
x = simulateData(current_chrom <- 6, growth_data, growth_traits, n_effects = 10)
sim_growth_data = x[[1]]
true_effects = x[[2]]

#growth_data[growth_traits] = scale(growth_data[growth_traits])
apply(as.matrix(sim_growth_data[growth_traits]), 2, var)

sim_model = runStanModel(current_chrom, growth_data, growth_traits,
                         chain = 3, iter = 350, warmup = 200, model_file = './SUR_horseShoePlus.stan')
chrom6_model = runStanModel(current_chrom, growth_data, growth_traits,
                         chain = 3, iter = 350, warmup = 200, model_file = './SUR_horseShoePlus.stan')

colMeans(sim_model[[3]]$w_ad)

png("./data/figures/lp.png")
plot(sim_model[[3]]$'tau_ad')
dev.off()

plotEffectEstimate(current_chrom, chrom6_model[[1]], "scaled_growth_0.1a_0.1ds_chrom6_HCp", scale = "free")
plotShrinkage(current_chrom, chrom6_model[[2]], "scaled_growth_0.1a_0.1ds_chrom6_HCp")

#plotEffectEstimate(current_chrom, sim_model[[1]], "scaled_growth_1s_chrom6_HCp", true_effects, scale = "free")
#plotShrinkage(current_chrom, sim_model[[2]], "scaled_growth_1s_chrom6_HCp")
