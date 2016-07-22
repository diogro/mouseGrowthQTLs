setwd("/home/diogro/projects/mouse-qtls")
source('stanFunctions.R')
source('read_mouse_data.R')
source('utils.R')

Rdatas_folder = "~/gdrive/LGSM_project_Rdatas/"
#Rdatas_folder = "./data/Rdatas/"

model_file = paste0(Rdatas_folder, "growth_scaled_HCPlus_stan")
readRDS()

c_model = chrom6_model
current_chrom = 6
treash = 0.4
shrink_test = colMeans(c_model[[3]]$shrink_ad) > treash | colMeans(c_model[[3]]$shrink_dm) > treash
selected_loci = which(apply(shrink_test, 2, any))

x = filter(c_model[[1]], marker %in% selected_loci)

plotEffectEstimate(6, x)

cols = c(growth_traits,
"FAMILY",
paste0("chrom", current_chrom, "_A", selected_loci),
paste0("chrom", current_chrom, "_D", selected_loci))

getStanInput(current_chrom,  growth_data[cols], growth_traits, length(selected_loci))

