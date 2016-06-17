setwd("/home/diogro/projects/mouse-qtls")
source('read_mouse_data.R')

makeSimData = function(data, marker, percent){
  a = findA(data, marker, percent)
  data + marker * a
}
markerVar = function(a, data, marker){
    data = mean(data) + residuals(lm(data ~ marker))
    za = data + marker * a
    n = length(data)
    V_pa = sum((za - mean(za))^2)/(n - 1)
    gen_means = tapply(za, marker, mean)
    gen_n = table(marker) / (n - 3)
    var_marker = gen_n %*% (gen_means - mean(za))^2
    (var_marker / V_pa)[1]
}
findA = function(data, marker, percent){
  V_p = var(data)
  a = sqrt((2*percent*V_p)/(1 - percent))
  return(a)
}

Rdatas_folder = "~/gdrive/LGSM_project_Rdatas/"
#Rdatas_folder = "./data/Rdatas/"

library(rstan)

area_data = inner_join(area_phen_std,
                       simulated_genomes[[2]],
                       by = "ID")

area_data = area_data %>% mutate_each(funs(scale), matches('area'))

true_effects = c(findA(area_data$area1, area_data$chrom4_A10, 0.05),
                 findA(area_data$area2, area_data$chrom4_A10, 0.05),
                 findA(area_data$area3, area_data$chrom4_A13, 0.05),
                 findA(area_data$area4, area_data$chrom4_A14, 0.05),
                 findA(area_data$area5, area_data$chrom4_A14, 0.05),
                 findA(area_data$area6, area_data$chrom4_A13, 0.05),
                 findA(area_data$area7, area_data$chrom4_A14, 0.05))
true_effects = data.frame(true_effects, trait = area_traits, type = "additive",
                          marker = c(rep(10, 2), 13, rep(14, 2), 13, 14))

area_data$area1 = makeSimData(area_data$area1, area_data$chrom4_A10, 0.05)
area_data$area2 = makeSimData(area_data$area2, area_data$chrom4_A10, 0.05)
area_data$area3 = makeSimData(area_data$area3, area_data$chrom4_A13, 0.05)
area_data$area4 = makeSimData(area_data$area4, area_data$chrom4_A14, 0.05)
area_data$area5 = makeSimData(area_data$area5, area_data$chrom4_A14, 0.05)
area_data$area6 = makeSimData(area_data$area6, area_data$chrom4_A13, 0.05)
area_data$area7 = makeSimData(area_data$area7, area_data$chrom4_A14, 0.05)

current_chrom = 4
getStanInput = function(current_chrom){
    K        = num_area_traits
    J        = loci_per_chrom[current_chrom]
    N        = dim(area_data)[1]
    n_family = length(unique(area_data$FAMILY))
    family   = as.integer(as.factor(area_data$FAMILY))
    ad       = as.matrix(select(area_data, matches('chrom4_A')))
    dm       = as.matrix(select(area_data, matches('chrom4_D')))
    y        = as.matrix(select(area_data, matches('area')))
    param_list = list(K        =  K,
                      J        =  J,
                      N        =  N,
                      n_family =  n_family,
                      family   =  family,
                      ad       =  ad,
                      dm       =  dm,
                      y        =  y)
    return(param_list)
}
stan_model = stan(file = './mixedModelGmatrix.stan',
                  data = getStanInput(4), chain=1, iter = 2000)
