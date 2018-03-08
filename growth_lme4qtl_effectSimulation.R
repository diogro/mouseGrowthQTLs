if(!require(purrr)){install.packages("purrr"); library(purrr)}
if(!require(corrplot)){install.packages("corrplot"); library(corrplot)}

setwd("/home/diogro/projects/mouse-qtls")
source('read_mouse_data.R')

#Rdatas_folder = "~/gdrive/LGSM_project_Rdatas/"
Rdatas_folder = "./Rdatas/"

vectorCor = function(x, y) Normalize(x) %*% Normalize(y)
markerCov = function(marker1, marker2){
  marker1_col = growth_markers[,makeMarkerList(marker1)]
  marker2_col = growth_markers[,makeMarkerList(marker2)]
  cov(marker1_col, marker2_col)[1]
}
makeMarkerList = function(pos) paste('chrom', pos[1],"_", 'A', pos[2], sep = '')
lt = function(x) x[lower.tri(x, diag = TRUE)]

install_load("doMC", "lme4qtl", "qvalue")
registerDoMC(15)

growth_data = inner_join(growth_phen_std,
                         growth_markers,
                         by = "ID") %>%
  gather(variable, value, growth12:growth78)
growth_data = mutate(growth_data, FAMILY = as.factor(FAMILY), variable = as.factor(variable))

load(file = "./Rdatas/significant_stan_fit.Rdata")
load(file = "Rdatas/growth_add_dom_effectsMatrix.Rdata")
significantMarkerMatrix = read_csv("./data/growth_significant_markers.csv")


simulateVg = function(a_effects, d_effects, markerMatrix){
  n_traits = dim(a_effects)[1]
  V_a = matrix(0, n_traits, n_traits)
  V_d = matrix(0, n_traits, n_traits)
  singleMakerSimulation = function(i){
    current_chrom = markerMatrix[i,1]
    focal_marker = markerMatrix[i,]
    # Variance due to focal marker
    p = rbeta(1, 1.1, 1.1)
    q = 1-p
    # additive contribution to Va
    V_a = 2*p*q * outer(a_effects[,i], a_effects[,i]) 
    # Variance due to LD with focal marker
    # diseqVa = function(j){
    #   if (i != j){ V_a = markerCov(focal_marker, markerMatrix[j,]) * outer(a_effects[,i], a_effects[,j])
    #   } else V_a = 0
    #   V_a
    # }
    # Va = V_a + Reduce("+", llply(1:nrow(markerMatrix), diseqVa))
    V_d = (2*p*q)^2 * outer(d_effects[,i], d_effects[,i]) 
    # Variance due to LD with focal marker
    # diseqVd = function(j){
    #   d2ij = markerCov(focal_marker, markerMatrix[j,])
    #   if (i != j){ V_d = d2ij^2 * outer(d_effects[,i], d_effects[,j])
    #   } else V_d = 0
    #   V_d
    # }
    # Vd = V_d + Reduce("+", llply(1:nrow(markerMatrix), diseqVd))
    V_a + V_d
  }
  Reduce("+", llply(1:dim(markerMatrix)[1], singleMakerSimulation, .parallel = TRUE))
}
simulateVg(effect_matrix_additive, effect_matrix_dominance, significantMarkerMatrix)
sim_Vg = rlply(100, simulateVg(effect_matrix_additive, effect_matrix_dominance, significantMarkerMatrix))
hist(RandomSkewers(sim_Vg, G)[,1])
