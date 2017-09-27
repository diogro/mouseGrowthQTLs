setwd("/home/diogro/projects/mouse-qtls")
source('stanFunctions.R')
source('read_mouse_data.R')
source('utils.R')

#Rdatas_folder = "~/gdrive/LGSM_project_Rdatas/"
Rdatas_folder = "./data/Rdatas/"

chroms = unlist(sapply(1:19, function(x) rep(x, loci_per_chrom[x])))
chrom_limits = data.frame(count = cumsum(loci_per_chrom) + 0.5)
chrom_limits$position = chrom_limits$count - c(31.5, diff(chrom_limits$count))/2
chrom_limits$chrom = 1:19 

model_file = paste0(Rdatas_folder, "necropsy_scaled_allmarkers_HCPlus")
necropsy_model = readRDS(model_file)
necropsy_shrink_out = rstan::summary( necropsy_model, pars = c("shrink_ad", "shrink_dm"))
necropsy_shrink = apply(matrix(necropsy_shrink_out[[1]][,'mean'], 10, 353, byrow = T), 2, max)

model_file = paste0(Rdatas_folder, "growth_scaled_allmarkers_HCPlus")
growth_model = readRDS(model_file)
growth_shrink_out = rstan::summary( growth_model, pars = c("shrink_ad", "shrink_dm"))
growth_shrink = apply(matrix(growth_shrink_out[[1]][,'mean'], 14, 353, byrow = T), 2, max)

model_file = paste0(Rdatas_folder, "area_scaled_allmarkers_HCPlus")
area_model = readRDS(model_file)
area_shrink_out = rstan::summary( area_model, pars = c("shrink_ad", "shrink_dm"))
area_shrink = apply(matrix(area_shrink_out[[1]][,'mean'], 14, 353, byrow = T), 2, max)

shrink_plots = list(
necropsy_shrink_plot =
    ggplot(data.frame(y = necropsy_shrink, 
                      x = seq(353),
                      chrom = chroms), aes(x, y)) + 
          geom_line() + 
          geom_hline(yintercept = 0.5, linetype = "dashed") + 
          geom_vline(data = chrom_limits, linetype = "dotted", aes(xintercept = count)) +
          labs(x = "marker", y = "Shrinkage intensity") +
          geom_text(data = chrom_limits, aes(x = position, y = 0.9, label = chrom)) + 
          ggtitle("Necropsy organ weights") + annotate("text", x = 0, y  = 1, label = "Chromossome"),
growth_shrink_plot =
    ggplot(data.frame(y = growth_shrink, 
                      x = seq(353),
                      chrom = chroms), aes(x, y)) + 
          geom_line() + 
          geom_hline(yintercept = 0.5, linetype = "dashed") + 
          geom_vline(data = chrom_limits, linetype = "dotted", aes(xintercept = count)) +
          ggtitle("Weekly growth") + labs(x = "marker", y = "Shrinkage intensity"),
area_shrink_plot =
    ggplot(data.frame(y = area_shrink, 
                      x = seq(353),
                      chrom = chroms), aes(x, y)) + 
          geom_line() + 
          geom_hline(yintercept = 0.5, linetype = "dashed") + 
          geom_vline(data = chrom_limits, linetype = "dotted", aes(xintercept = count)) +
          ggtitle("Mandible areas") + labs(x = "marker", y = "Shrinkage intensity")
)      


save_plot("~/images/hs_shrink.png", shrink_plot, base_height = 5, base_aspect_ratio = 2)

shrink_plot = plot_grid(shrink_plots[[1]], 
                        shrink_plots[[2]], 
                        shrink_plots[[3]], ncol = 1)
save_plot("~/images/hs_shrink.png", shrink_plot, nrow = 3, base_height = 2.5, base_aspect_ratio = 5)

shrink_list = list(necropsy = necropsy_shrink_out[[1]], growth = growth_shrink_out[[1]], area = area_shrink_out[[1]])
shrink_list = list(necropsy = necropsy_shrink, growth = growth_shrink, area = area_shrink)
saveRDS(shrink_list, file = paste0(Rdatas_folder, "shrink_mapping"))

# Genetic effects
chroms = unlist(sapply(1:19, function(x) rep(x, loci_per_chrom[x])))
markers_vector =  unlist(sapply(1:19, function(x) seq(loci_per_chrom[x])))

model_file = paste0(Rdatas_folder, "necropsy_scaled_allmarkers_HCPlus")
necropsy_model = readRDS(model_file)
necropsy_effects = getStanEffects(chroms, necropsy_model, necropsy_traits, 353, markers_vector)

model_file = paste0(Rdatas_folder, "growth_scaled_allmarkers_HCPlus")
growth_model = readRDS(model_file)
growth_effects = getStanEffects(chroms, growth_model, growth_traits, 353, markers_vector)

model_file = paste0(Rdatas_folder, "area_scaled_allmarkers_HCPlus")
area_model = readRDS(model_file)
area_effects = getStanEffects(chroms, area_model, area_traits, 353, markers_vector)

effect_list = list(necropsy = necropsy_effects, growth = growth_effects, area = area_effects)
saveRDS(effect_list, file = paste0(Rdatas_folder, "hsplus_effects"))
