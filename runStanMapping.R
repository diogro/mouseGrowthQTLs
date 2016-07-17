setwd("/home/diogro/projects/mouse-qtls")
source("stanFunctions.R")
source('read_mouse_data.R')
source('utils.R')

#Rdatas_folder = "~/gdrive/LGSM_project_Rdatas/"
Rdatas_folder = "./data/Rdatas/"

weight_data = inner_join(weight_phen_std, weight_markers, by = "ID")

weight_data[weight_traits] = scale(weight_data[weight_traits])

#current_data = weight_data
#trait_vector = weight_traits
#current_chrom = 19
#stan_model = stan_model_SUR_HC
#iter = 10
#chain = 1
#parallel = TRUE

library(doMC)
registerDoMC(19)
weight_mapping = runStanModelFullGenome(weight_data, weight_traits, 2000, TRUE)
ef_plots = llply(1:19, plotEffectEstimate, weight_mapping[[1]], "weight")
ef_plots = llply(1:19, plotWeights, weight_mapping[[2]], "weight")

