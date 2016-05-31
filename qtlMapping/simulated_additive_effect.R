setwd("~/projects/mouse-qtls/")
source("~/projects/mouse-qtls/read_mouse_data.R")
source("~/projects/mouse-qtls/OAuth_lem_server.R")

a = 2.2
markerVar = function(x){
    data$za = data$z + data$marker * a/2
    V_pa = sum((data$za - mean(data$za))^2)/(1000 - 1)
    gen_means = tapply(data$za, data$marker, mean)
    gen_n = table(data$marker)
    var_marker = gen_n %*% (gen_means - mean(data$za))^2 / (1000 - 3)
    (var_marker / V_pa)[1]
}

simulated_markers


