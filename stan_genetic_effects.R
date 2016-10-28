setwd("/home/diogro/projects/mouse-qtls")
source('stanFunctions.R')
source('read_mouse_data.R')
source('utils.R')

#Rdatas_folder = "~/gdrive/LGSM_project_Rdatas/"
Rdatas_folder = "./data/Rdatas/"

model_file = paste0(Rdatas_folder, "growth_scaled_HCPlus_stan")
c_model = readRDS(model_file)

chroms = unlist(sapply(1:19, function(x) rep(x, loci_per_chrom[x])))

shrink_out = rstan::summary( c_model, pars = c("shrink_ad", "shrink_dm"))
ggplot(data.frame(y = as.numeric(c_model[[2]][,"mean"][[1]]), 
                  x = seq(353),
                  trait = rep(growth_traits, each = 353),
                  chrom = chroms), aes(x, y)) + 
  geom_line() + geom_point(size = 1) + 
  geom_hline(yintercept = 0.5, linetype = "dashed") + facet_wrap(~trait, ncol = 1) 

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

