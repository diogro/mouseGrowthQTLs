setwd("/home/diogro/projects/mouse-qtls")
source("stanFunctions.R")
source('read_mouse_data.R')
source('utils.R')

#Rdatas_folder = "~/gdrive/LGSM_project_Rdatas/"
Rdatas_folder = "./data/Rdatas/"

#current_data = weight_data
#trait_vector = weight_traits
#current_chrom = 19
#stan_model = stan_model_SUR_HC
#iter = 10
#chain = 1
#parallel = TRUE

install_load("doMC")
registerDoMC(19)

growth_data = inner_join(growth_phen_std, growth_markers, by = "ID")
growth_data[growth_traits] = scale(growth_data[growth_traits])
growth_mapping = runStanModelFullGenome(growth_data, growth_traits, 500, 100, TRUE)
ef_plots = llply(1:19, plotEffectEstimate, growth_mapping[[1]], "growth")
ef_plots = llply(1:19, plotShrinkage, growth_mapping[[2]], "growth")
model_file = paste0(Rdatas_folder, "growth_scaled_HCPlus_stan")
saveRDS(growth_mapping, model_file)

area_data = inner_join(area_phen_std, area_markers, by = "ID")
area_data[area_traits] = scale(area_data[area_traits])
area_mapping = runStanModelFullGenome(area_data, area_traits, 500, 100, TRUE)
ef_plots = llply(1:19, plotEffectEstimate, area_mapping[[1]], "area")
ef_plots = llply(1:19, plotShrinkage, area_mapping[[2]], "area")
model_file = paste0(Rdatas_folder, "area_scaled_HCPlus_stan")
saveRDS(area_mapping, model_file)

necropsy_data = inner_join(necropsy_phen_std, necropsy_markers, by = "ID")
necropsy_data[necropsy_traits] = scale(necropsy_data[necropsy_traits])
necropsy_mapping = runStanModelFullGenome(necropsy_data, necropsy_traits, 500, 100, TRUE)
ef_plots = llply(1:19, plotEffectEstimate, necropsy_mapping[[1]], "necropsy")
ef_plots = llply(1:19, plotShrinkage, necropsy_mapping[[2]], "necropsy")
model_file = paste0(Rdatas_folder, "necropsy_scaled_HCPlus_stan")
saveRDS(necropsy_mapping, model_file)
