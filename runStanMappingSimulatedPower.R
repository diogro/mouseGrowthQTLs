setwd("/home/diogro/projects/mouse-qtls")
source('stanFunctions.R')
source('read_mouse_data.R')
source('utils.R')

#Rdatas_folder = "~/gdrive/LGSM_project_Rdatas/"
Rdatas_folder = "./data/Rdatas/"

options(mc.cores = 4)

current_genome = sample(length(simulated_genomes), 1)
current_effect = sample(1:353, 1)
focal = as.matrix(simulated_genomes[[current_genome]] %>% select(matches("_A")))[,current_effect]
cor_dist = lapply(simulated_markers,
                  function(x) apply(as.matrix(x %>% select(matches("_A"))),
                                    2, function(x) cor(x, focal)))
max(unlist(cor_dist))
sort(cor_dist, d = T)[2]
quantile(cor_dist, 0.999999)

lx = x[lower.tri(x)]
summary(lx)

area_data = inner_join(area_phen_std,
                       simulated_genomes[[5]],
                       by = "ID")

area_data[area_traits] = scale(area_data[area_traits])
area_data[area_traits] = scale(area_data[area_traits], scale = rep(sqrt(1/40), num_area_traits))

true_effects = c(findA(area_data$area1, area_data$chrom4_A10, 0.02),
                 findA(area_data$area2, area_data$chrom4_A10, 0.02),
                 findA(area_data$area3, area_data$chrom4_A13, 0.02),
                 findA(area_data$area4, area_data$chrom4_A14, 0.02),
                 findA(area_data$area5, area_data$chrom4_A14, 0.02),
                 findA(area_data$area6, area_data$chrom4_A13, 0.02),
                 findA(area_data$area7, area_data$chrom4_A14, 0.02))
true_effects = data.frame(true_effects, trait = area_traits, type = "aditive",
                          marker = c(rep(10, 2), 13, rep(14, 2), 13, 14), chrom = 4)

area_data$area1 = makeSimData(area_data$area1, area_data$chrom4_A10, 0.02)
area_data$area2 = makeSimData(area_data$area2, area_data$chrom4_A10, 0.02)
area_data$area3 = makeSimData(area_data$area3, area_data$chrom4_A13, 0.02)
area_data$area4 = makeSimData(area_data$area4, area_data$chrom4_A14, 0.02)
area_data$area5 = makeSimData(area_data$area5, area_data$chrom4_A14, 0.02)
area_data$area6 = makeSimData(area_data$area6, area_data$chrom4_A13, 0.02)
area_data$area7 = makeSimData(area_data$area7, area_data$chrom4_A14, 0.02)

area_data[area_traits] = scale(area_data[area_traits])
area_data[area_traits] = scale(area_data[area_traits], scale = rep(sqrt(1/40), num_area_traits))
apply(as.matrix(area_data[area_traits]), 2, var)

current_chrom = 4
sim_model = runStanModel(current_chrom, area_data, area_traits, chain = 4, iter = 400, model_file = './SUR_horseShoePlus.stan', 
                         control = list(adapt_delta = 0.95))

plotEffectEstimate(current_chrom, sim_model[[1]], "simPlus_area_scaled_2pct", true_effects)
plotWeights(current_chrom, sim_model[[2]], "simPlus_area_scaled_2pct")
