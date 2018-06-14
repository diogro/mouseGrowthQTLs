setwd("/home/diogro/projects/mouseGrowthQTLs/")
source("stanFunctions.R")
source('read_mouse_data.R')
source('utils.R')

#Rdatas_folder = "~/gdrive/LGSM_project_Rdatas/"
Rdatas_folder = "./data/Rdatas/"

rstan_options(auto_write = TRUE)
options(mc.cores = 4)

growth_data = inner_join(growth_phen_std, growth_markers, by = "ID")
growth_data[growth_traits] = scale(growth_data[growth_traits])
growth_full_genome = stan(file = './SUR_horseShoe.stan',
                          data = getStanInput(1, growth_data, growth_traits),
                          chain=1, iter = 200, 
                          control = list(adapt_delta=0.9, max_treedepth=15))

#model_file = paste0(Rdatas_folder, "growth_scaled_allmarkers_HCPlus")
#saveRDS(growth_full_genome, model_file)
betas = getStanEffects(1, growth_full_genome, growth_traits)
plotEffectEstimate(1, betas)

