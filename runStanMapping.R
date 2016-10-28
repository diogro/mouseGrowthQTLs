setwd("/home/diogro/projects/mouse-qtls")
source("stanFunctions.R")
source('read_mouse_data.R')
source('utils.R')

#Rdatas_folder = "~/gdrive/LGSM_project_Rdatas/"
Rdatas_folder = "./data/Rdatas/"

rstan_options(auto_write = TRUE)
options(mc.cores = 4)

growth_data = inner_join(growth_phen_std, growth_markers, by = "ID")
growth_data[growth_traits] = scale(growth_data[growth_traits])
growth_full_genome = stan(file = './SUR_horseShoePlus.stan',
                          data = getStanInputFullGenome(growth_data, growth_traits),
                          chain=4, iter = 6, warmup = 2)
model_file = paste0(Rdatas_folder, "growth_scaled_allmarkers_HCPlus")
saveRDS(growth_full_genome, model_file)

area_data = inner_join(area_phen_std, area_markers, by = "ID")
area_data[area_traits] = scale(area_data[area_traits])
area_full_genome = stan(file = './SUR_horseShoePlus.stan',
                          data = getStanInputFullGenome(current_data = area_data, trait_vector = area_traits),
                          chain=4, iter = 600, warmup = 200)
model_file = paste0(Rdatas_folder, "area_scaled_allmarkers_HCPlus")
saveRDS(area_full_genome, model_file)

necropsy_data = inner_join(necropsy_phen_std, necropsy_markers, by = "ID")
necropsy_data[necropsy_traits] = scale(necropsy_data[necropsy_traits])
necropsy_full_genome = stan(file = './SUR_horseShoePlus.stan',
                          data = getStanInputFullGenome(current_data = necropsy_data, trait_vector = necropsy_traits),
                          chain=4, iter = 600, warmup = 200)
model_file = paste0(Rdatas_folder, "necropsy_scaled_allmarkers_HCPlus")
saveRDS(necropsy_full_genome, model_file)