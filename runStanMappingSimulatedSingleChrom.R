setwd("/home/diogro/projects/mouse-qtls")
source('stanFunctions.R')
source('read_mouse_data.R')
source('utils.R')

#Rdatas_folder = "~/gdrive/LGSM_project_Rdatas/"
Rdatas_folder = "./data/Rdatas/"

options(mc.cores = 4)

weight_data = inner_join(weight_phen_std,
                       simulated_genomes[[8]],
                       by = "ID")

weight_data = inner_join(weight_phen_std,
                       markers,
                       by = "ID")

weight_data[weight_traits] = scale(weight_data[weight_traits])
weight_data[weight_traits] = scale(weight_data[weight_traits], scale = rep(sqrt(1/20), num_weight_traits))

true_effects = c(findA(weight_data$WEEK1, weight_data$chrom4_A1, 0.02),
                 findA(weight_data$WEEK2, weight_data$chrom4_A10, 0.02),
                 findA(weight_data$WEEK3, weight_data$chrom4_A13, 0.02),
                 findA(weight_data$WEEK4, weight_data$chrom4_A14, 0.02),
                 findA(weight_data$WEEK5, weight_data$chrom4_A14, 0.02),
                 findA(weight_data$WEEK6, weight_data$chrom4_A13, 0.02),
                 findA(weight_data$WEEK7, weight_data$chrom4_A14, 0.02))
true_effects = data.frame(true_effects, trait = weight_traits[1:7], type = "aditive",
                          marker = c(1, 10, 13, rep(14, 2), 13, 14), chrom = 4)

weight_data$WEEK1 = makeSimData(weight_data$WEEK1, weight_data$chrom4_A1, 0.02)
weight_data$WEEK2 = makeSimData(weight_data$WEEK2, weight_data$chrom4_A10, 0.02)
weight_data$WEEK3 = makeSimData(weight_data$WEEK3, weight_data$chrom4_A13, 0.02)
weight_data$WEEK4 = makeSimData(weight_data$WEEK4, weight_data$chrom4_A14, 0.02)
weight_data$WEEK5 = makeSimData(weight_data$WEEK5, weight_data$chrom4_A14, 0.02)
weight_data$WEEK6 = makeSimData(weight_data$WEEK6, weight_data$chrom4_A13, 0.02)
weight_data$WEEK7 = makeSimData(weight_data$WEEK7, weight_data$chrom4_A14, 0.02)

weight_data[weight_traits] = scale(weight_data[weight_traits])
#weight_data[weight_traits] = scale(weight_data[weight_traits], scale = rep(sqrt(1/40), num_weight_traits))
apply(as.matrix(weight_data[weight_traits]), 2, var)

current_chrom = 6
sim_model = runStanModel(current_chrom, weight_data, weight_traits,
                         chain = 4, iter = 200, model_file = './SUR_horseShoePlus.stan',
                         control = list(adapt_delta = 0.95))

plotEffectEstimate(current_chrom, sim_model[[1]], "weight", true_effects)
plotWeights(current_chrom, sim_model[[2]], "weight")

plotEffectEstimate(current_chrom, sim_model[[1]], "weight_0.2s_2pct", true_effects)
plotWeights(current_chrom, sim_model[[2]], "simPlus_weight_0.2s_2pct")
