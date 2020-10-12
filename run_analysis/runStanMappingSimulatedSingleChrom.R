setwd("/home/diogro/projects/mouse-qtls")
source('stanFunctions.R')
source('read_mouse_data.R')
source('utils.R')

#Rdatas_folder = "~/gdrive/LGSM_project_Rdatas/"
Rdatas_folder = "./data/Rdatas/"

rstan_options(auto_write = TRUE)
options(mc.cores = 3)

growth_data = inner_join(growth_phen_std, simulated_genomes[[20]], by = "ID")
growth_data[growth_traits] = scale(growth_data[growth_traits])

#growth_data = inner_join(growth_phen_std, markers, by = "ID")
set.seed(39)
simulateData = function(current_chrom, current_data, trait_vector,
                        n_effects = length(trait_vector), eff_sd = 0.5){
    true_effects = NULL
    all_loci = sample(loci_per_chrom[current_chrom], n_effects)
    for(i in 1:n_effects){
        current_trait = sample(trait_vector, sample(1:4, 1))
        type = sample(c("A", "D"), 1)
        current_loci = all_loci[i]
        locus = paste0("chrom", current_chrom, "_", type, current_loci)
        effect = 0.3#rnorm(1, 0, eff_sd)
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

true_effects$plot_position =  31 * 7 * 2 - (31 * (as.numeric(factor(true_effects$trait, levels = growth_traits)) - 1) + true_effects$marker + 31 * 7 * as.numeric(true_effects$type == "dominance"))

#growth_data[growth_traits] = scale(growth_data[growth_traits])
apply(as.matrix(sim_growth_data[growth_traits]), 2, var)

# stan_model_SUR_HC = stan(file = './SUR_horseShoe.stan',
#                          data = getStanInput(current_chrom, sim_growth_data, growth_traits),
#                          chain=3, iter = 200, warmup = 100)
stan_input = getStanInput(current_chrom, sim_growth_data, growth_traits)
stan_input$ad = stan_input$ad / 5
stan_input$dm = stan_input$dm / 5
stan_model_SUR_HCp = stan(file = './SUR_horseShoePlus.stan',
                          data = stan_input,
                          chain=3, iter = 200, warmup = 100)
# stan_model_SUR_HAL = stan(file = './SUR_HAL.stan',
#                           data = getStanInput(current_chrom, sim_growth_data, growth_traits),
#                           chain=3, iter = 200, warmup = 100)
# stan_model_SUR = stan(file = './SUR.stan',
#                       data = getStanInput(current_chrom, sim_growth_data, growth_traits),
#                       chain=3, iter = 200, warmup = 100)
effects = getStanEffects(current_chrom, stan_model_SUR_HCp, growth_traits)
weights = getStanShrinkage(current_chrom, stan_model_SUR_HCp, growth_traits)
sim_model = (list(effects, weights, rstan::extract(stan_model_SUR_HCp)))

true_effects$true_effects = true_effects$true_effects * 5
plotEffectEstimate(current_chrom, effects, ,true_effects)
plotShrinkage(current_chrom, weights)

# x = plot( stan_model_SUR_HC, pars = c("w_ad", "w_dm")) + theme(axis.text.y=element_blank()) +
# geom_vline(xintercept = 0.2, linetype = "dashed")+
# geom_vline(xintercept = 0)+
# scale_x_continuous(lim = c(-0.7, 0.7))
z = plot( stan_model_SUR_HCp, pars = c("w_ad", "w_dm")) + theme(axis.text.y=element_blank(),
                                                                axis.title = element_text(size = 14)) +
geom_vline(xintercept = 0.2, linetype = "dashed", size = 1) +
geom_vline(xintercept = 0) + geom_hline(data = true_effects, aes(yintercept = plot_position), color = "blue") +
labs(x = "Effect size", y = "Coefficient index")
# w = plot( stan_model_SUR_HAL, pars = c("w_ad", "w_dm")) + theme(axis.text.y=element_blank()) +
# geom_vline(xintercept = 0.2, linetype = "dashed")+
# geom_vline(xintercept = 0)+
# scale_x_continuous(lim = c(-0.7, 0.7))
y = plot( stan_model_SUR, pars = c("w_ad", "w_dm")) + theme(axis.text.y=element_blank(),
                                                            axis.title = element_text(size = 14)) +
geom_vline(xintercept = 0.2, linetype = "dashed", size = 1) +
geom_vline(xintercept = 0) + geom_hline(data = true_effects, aes(yintercept = plot_position), color = "blue") +
labs(x = "Effect size", y = "Coefficient index")
# plot_grid(y, w, x, z, ncol = 4, labels = c("ML", "HAL", "Horseshoe", "Horseshoe+"))
save_plot("~/images/normal_prior.png", y, base_height = 5)
yz = plot_grid( y, z, ncol = 2, labels = c("Guassian", "Horseshoe+"))
save_plot("~/images/hs_prior.png", yz, base_height = 5, ncol = 2)

shrink_out = rstan::summary( stan_model_SUR_HCp, pars = c("shrink_ad", "shrink_dm"))
shrink_plot = ggplot(data.frame(y = shrink_out[[1]][,"mean"], x = 1:length(shrink_out[[1]][,"mean"])), aes(x, y)) + 
  geom_line() + geom_point(size = 2) + geom_hline(yintercept = 0.5, linetype = "dashed") + 
  geom_vline(data = true_effects, aes(xintercept = 7 * 31 * 2 - plot_position), color = "blue", linetype = "dotted") +
  labs(x = "Coefficient index", y = "Shrinkage")

save_plot("~/images/sim_shrinkage.png", shrink_plot, base_height = 5, base_aspect_ratio = 1.8)
