setwd("/home/diogro/projects/mouse-qtls")
source("stanFunction.R")
source('read_mouse_data.R')
source('utils.R')

#Rdatas_folder = "~/gdrive/LGSM_project_Rdatas/"
Rdatas_folder = "./data/Rdatas/"

area_data = inner_join(area_phen_std,
                       Reduce(inner_join, area_markers),
                       by = "ID")

area_data = area_data %>% mutate_each(funs(scale), matches('area'))
area_data[area_traits] = scale(area_data[area_traits])
area_data[area_traits] = scale(area_data[area_traits], scale = rep(sqrt(1/40), num_area_traits))

#current_data = area_data
#trait_vector = area_traits
#current_chrom = 19
#stan_model = stan_model_SUR_HC
#iter = 10
#chain = 1
#parallel = TRUE

library(doMC)
registerDoMC(19)
area_mapping = runStanModelFullGenome(area_data, area_traits, 2000, TRUE)
ef_plots = llply(1:19, plotEffectEstimate, area_mapping[[1]], "area_scaled")
ef_plots = llply(1:19, plotWeights, area_mapping[[2]], "area_scaled")
