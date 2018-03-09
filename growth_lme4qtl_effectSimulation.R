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
  trait_sd = sapply(growth_phen[,growth_traits], sd)
  a_effects = a_effects * trait_sd
  d_effects = d_effects * trait_sd
  singleMakerSimulation = function(i){
    current_chrom = markerMatrix[i,1]
    focal_marker = markerMatrix[i,]
    # Variance due to focal marker
    p = rbeta(1, 1.1, 1.1)
    q = 1-p
    # additive contribution to Va
    V_a = 2*p*q * outer(a_effects[,i], a_effects[,i]) 
    # Dominance contribution to Va
    V_a = V_a + 2*p*q * (q - p)^2 * outer(d_effects[,i], d_effects[,i])
    # Additive by dominance contribution to Va
    V_a = V_a + 2*p*q * (q - p)  * (outer(a_effects[,i], d_effects[,i]) + 
                                    outer(d_effects[,i], a_effects[,i]))
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
G
sim_Vg = rlply(100, simulateVg(effect_matrix_additive, effect_matrix_dominance, significantMarkerMatrix))

Comparisons = data.frame(Krzanowski = KrzCor(sim_Vg, G),
           RandomSkewers = RandomSkewers(sim_Vg, G)[,1]) %>%  gather() 
  
ggplot(Comparisons, aes(value, group = key, fill = key)) + geom_density(alpha = 0.5) + 
  scale_x_continuous(limits = c(0.7, 1.0)) + 
  geom_vline(xintercept = RandomSkewers(Vg_mean, G)[1], color = "seagreen") + 
  geom_vline(xintercept = KrzCor(Vg_mean, G), color = "red") + 
  labs(x = "Matrix similarity with Full-Sib Family G Matrix", y = "Density") + 
  annotate("text", x = 0.955, y = 22, label = "Comparison with observed\n allele frequencies")

ggplot(ldply(sim_Vg, CalcInt), aes(V1)) + geom_density(fill = "blue", alpha = 0.5) + 
  geom_vline(xintercept = CalcInt(Vg_mean), color = "red")

CalcInt(G)
CalcInt(Vg_mean)
corrG = cov2cor(G)
corrVg = cov2cor(Vg_mean)
plot(lt(G, FALSE) ~ lt(Vg_mean, FALSE))
abline(0, 1)
