setwd("/home/diogro/projects/mouse-qtls")
source('stanFunctions.R')
source('read_mouse_data.R')
source('utils.R')

Rdatas_folder = "~/gdrive/LGSM_project_Rdatas/"
#Rdatas_folder = "./data/Rdatas/"

rstan_options(auto_write = TRUE)
options(mc.cores = 3)

growth_data = inner_join(growth_phen_std, simulated_genomes[[20]], by = "ID")
growth_data[growth_traits] = scale(growth_data[growth_traits])

#growth_data = inner_join(growth_phen_std, markers, by = "ID")
simulateData = function(current_chrom, current_data, trait_vector,
                        n_effects = length(trait_vector), eff_sd = 0.5){
    true_effects = NULL
    all_loci = sample(loci_per_chrom[current_chrom], n_effects)
    for(i in 1:n_effects){
        current_trait = sample(trait_vector, sample(1:4, 1))
        type = sample(c("A", "D"), 1)
        current_loci = all_loci[i]
        locus = paste0("chrom", current_chrom, "_", type, current_loci)
        effect = 0.2#rnorm(1, 0, eff_sd)
        for(trait in current_trait){
            current_data[trait] = current_data[trait] + effect * current_data[locus]
            true_effects = rbind(true_effects, data.frame(chrom = current_chrom, marker = current_loci,
                                                          type = ifelse(type == "A", "additive", "dominance"),
                                                          true_effects = effect, trait = trait))
        }
    }
    return(list(current_data, true_effects))
}
x = simulateData(current_chrom <- 1, growth_data, growth_traits, n_effects = 3, eff_sd = 0.5)
sim_growth_data = x[[1]]
true_effects = x[[2]]

#growth_data[growth_traits] = scale(growth_data[growth_traits])
apply(as.matrix(sim_growth_data[growth_traits]), 2, var)

stan_model_SUR_HC = stan(file = './SUR_horseShoe.stan',
                         data = getStanInput(current_chrom, sim_growth_data, growth_traits),
                         chain=3, iter = 200, warmup = 100)
stan_model_SUR_HCp = stan(file = './SUR_horseShoePlus.stan',
                         data = getStanInput(current_chrom, sim_growth_data, growth_traits),
                         chain=3, iter = 200, warmup = 100)
stan_model_SUR_HAL = stan(file = './SUR_HAL.stan',
                         data = getStanInput(current_chrom, sim_growth_data, growth_traits),
                         chain=3, iter = 200, warmup = 100)
stan_model_SUR = stan(file = './SUR.stan',
                         data = getStanInput(current_chrom, sim_growth_data, growth_traits),
                         chain=3, iter = 200, warmup = 100)
effects = getStanEffects(current_chrom, stan_model_SUR_HAL, growth_traits)
weights = getStanShrinkage(current_chrom, stan_model_SUR_HAL, growth_traits)
sim_model = (list(effects, weights, rstan::extract(stan_model_SUR_HC)))

plotEffectEstimate(current_chrom, effects, ,true_effects)
plotShrinkage(current_chrom, weights)

x = plot( stan_model_SUR_HC, pars = c("w_ad", "w_dm")) + theme(axis.text.y=element_blank()) +
geom_vline(xintercept = 00.15, linetype = "dashed")+
geom_vline(xintercept = 0)+
scale_x_continuous(lim = c(-0.7, 0.7))
z = plot( stan_model_SUR_HCp, pars = c("w_ad", "w_dm")) + theme(axis.text.y=element_blank()) +
geom_vline(xintercept = 00.15, linetype = "dashed")+
geom_vline(xintercept = 0)+
scale_x_continuous(lim = c(-0.7, 0.7))
w = plot( stan_model_SUR_HAL, pars = c("w_ad", "w_dm")) + theme(axis.text.y=element_blank()) +
geom_vline(xintercept = 00.15, linetype = "dashed")+
geom_vline(xintercept = 0)+
scale_x_continuous(lim = c(-0.7, 0.7))
y = plot( stan_model_SUR, pars = c("w_ad", "w_dm")) + theme(axis.text.y=element_blank())+
geom_vline(xintercept = 00.15, linetype = "dashed")+
geom_vline(xintercept = 0)+
scale_x_continuous(lim = c(-0.7, 0.7))
plot_grid(y, w, x, z, ncol = 4, labels = c("ML", "HAL", "Horseshoe", "Horseshoe+"))
plot_grid( w, z, ncol = 2, labels = c("HAL", "Horseshoe+"))

plot( stan_model_SUR_HAL, pars = c("lambda_ad"))
plot( stan_model_SUR_HC, pars = c("lambda_ad"))
plot( stan_model_SUR_HCp, pars = c("etaLambda_ad"))
plot( stan_model_SUR_HAL, pars = c("shrink_ad"))
plot( stan_model_SUR_HC, pars = c("shrink_ad"))
plot( stan_model_SUR_HCp, pars = c("shrink_ad"))

traceplot( stan_model_SUR_HC,
  pars = c("shrink_dm"),
  inc_warmup = TRUE
)

sim_model = runStanModel(current_chrom, sim_growth_data, growth_traits,
                         chain = 3, iter = 450, warmup = 300, model_file = './SUR_horseShoePlus.stan')
chrom6_model = runStanModel(current_chrom, growth_data, growth_traits,
                         chain = 3, iter = 350, warmup = 200, model_file = './SUR_horseShoePlus.stan')

colMeans(sim_model[[3]]$w_ad)

png("./data/figures/lp.png")
plot(sim_model[[3]]$'tau_ad')
dev.off()

#plotEffectEstimate(current_chrom, chrom6_model[[1]], "scaled_growth_0.1a_0.1ds_chrom6_HCp", scale = "free")
#plotShrinkage(current_chrom, chrom6_model[[2]], "scaled_growth_0.1a_0.1ds_chrom6_HCp")

plotEffectEstimate(current_chrom, sim_model[[1]], "ss_growth_1s_chrom6_HCp", true_effects)
plotShrinkage(current_chrom, sim_model[[2]], "ss_growth_1s_chrom6_HCp")
