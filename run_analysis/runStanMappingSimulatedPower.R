setwd("/home/diogro/projects/mouse-qtls")
source('stanFunctions.R')
source('read_mouse_data.R')
source('utils.R')

Rdatas_folder = "~/gdrive/LGSM_project_Rdatas/"
#Rdatas_folder = "./data/Rdatas/"

rstan_options(auto_write = TRUE)
options(mc.cores = 3)

current_marker = 1
true_effect = 0.1
trait_vector = growth_traits
current_std_data = growth_phen_std

simulateDataSingleMarker = function(current_marker, true_effect, current_std_data, trait_vector){
    current_data = inner_join(current_std_data,
                             simulated_markers[[1]][c("ID",
                                                      paste0("chrom1_A", current_marker:10),
                                                      paste0("chrom1_D", current_marker:10))],
                             by = "ID")
    current_data[trait_vector] = scale(current_data[trait_vector])
    for(trait in trait_vector){
        current_data[trait] = current_data[trait] +
                              true_effect * current_data[paste0("chrom1_A", current_marker)] +
                              true_effect * current_data[paste0("chrom1_D", current_marker)]
    }
    return(current_data)
}
sim_data = simulateDataSingleMarker(1, 0.0, growth_phen_std, growth_traits)

stan_model_SUR_HC = stan(file = './SUR_horseShoe.stan',
                         data = getStanInput(1, sim_data, growth_traits, 10),
                         chain=3, iter = 200, warmup = 100,
                         control = list(adapt_delta = 0.95, stepsize = 0.01))
effects = getStanEffects(1, stan_model_SUR_HC, growth_traits, J  = 10, markers = 1:10)
weights = getStanShrinkage(1, stan_model_SUR_HC, growth_traits, J = 10, markers = 1:10)
weights %>% print(n = 20)
sim_model = (list(effects, weights, rstan::extract(stan_model_SUR_HC)))

stan_model_SUR = stan(file = './SUR.stan',
                         data = getStanInput(1, sim_data, growth_traits, 10),
                         chain=3, iter = 200, warmup = 100,
                         control = list(adapt_delta = 0.95, stepsize = 0.01))
effects_SUR = getStanEffects(1, stan_model_SUR, growth_traits, J  = 10, markers = 1:10)


traceplot( stan_model_SUR_HC, pars = c("shrink_dm"), inc_warmup = TRUE)
traceplot( stan_model_SUR_HC,
  pars = c("w_dm[1,1]","w_dm[1,10]", "w_dm[3,11]"),
  inc_warmup = F
)
traceplot( stan_model_SUR_HC, pars = c("tau_ad"), inc_warmup = TRUE)
traceplot( stan_model_SUR_HC, pars = c("w_ad"), inc_warmup = F)
traceplot( stan_model_SUR_HC, pars = c("w_ad[7, 14]"), inc_warmup = F)
plot( stan_model_SUR_HC, pars = c("shrink_ad"))
<cr>x = plot( stan_model_SUR_HC, pars = c("w_ad", "w_dm"))
y = plot( stan_model_SUR, pars = c("w_ad", "w_dm"))
plot_grid(x, y)
plot( stan_model_SUR_HC, pars = c("w_dm"))
plot( stan_model_SUR_HC, pars = c("tau_ad", "tau_dm"))
plot( stan_model_SUR_HC, pars = c("lambda_ad", "lambda_dm"))
plotEffectEstimate(1, effects)
plotEffectEstimate(1, effects_SUR)
plotShrinkage(1, weights)
