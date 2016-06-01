setwd("~/projects/mouse-qtls/")
source("~/projects/mouse-qtls/read_mouse_data.R")
source("~/projects/mouse-qtls/OAuth_lem_server.R")

data = data.frame(z = rnorm(1000, 1, 10), 
                  marker = sample(c(-1, 0, 1), 1000, replace = TRUE))

a = 2.2
markerVar = function(x){
    data$za = data$z + data$marker * a/2
    V_pa = sum((data$za - mean(data$za))^2)/(1000 - 1)
    gen_means = tapply(data$za, data$marker, mean)
    gen_n = table(data$marker) / (1000 - 3)
    var_marker = gen_n %*% (gen_means - mean(data$za))^2
    (var_marker / V_pa)[1]
}
markerVar(1)

V_p = sum((data$za - mean(data$za))^2)/(1000 - 1)
gen_means = tapply(data$za, data$marker, mean)
gen_n = table(data$marker)/1000
V_a = gen_n %*% (gen_means - mean(data$za))^2
(var_marker / V_pa)[1]

eff_var = 0.1

2*(eff_var * V_p / (1-eff_var) - V_a)/(gen_n[3] - gen_n[1])

sim_chrom_number = 1
area_data = inner_join(area_phen_std, simulated_markers[[sim_chrom_number]], by = "ID")

current_loci = 1
names(area_data)

sim_data = dplyr::select(area_data, area1:area7, matches("A1$"))

