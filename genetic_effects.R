source('./read_mouse_data.R')
if(!require(install.load)) {install.packages("install.load"); library(install.load)}
install_load("plyr", "evolqg", "dplyr", "tidyr", "readr", "ggplot2", "cowplot")

nullDIC    = readRDS(file = paste0("./data/area\ traits/nullDIC.rds"))
nullDIC_df = reshape2::melt(nullDIC)
names(nullDIC_df)
quantile(nullDIC_df$value, 0.999)
nullDIC_plot = ggplot(nullDIC_df, aes(value)) + geom_histogram(binwidth = 0.3) + geom_vline(xintercept = quantile(nullDIC_df$value, 0.999)) + annotate("text", x = 16, y = 300, label = "99.9 percentile", hjust = 0) + labs(x = "DIC difference") + scale_x_continuous(breaks = c(-20, -10, 0, 10, 15, 20))
save_plot("~/Dropbox/labbio/relatorios/fapesp/fapesp-relatorio-2016-10-30-BEPE/images/null_DIC.png", nullDIC_plot, base_height = 4, base_aspect_ratio = 1.8)


getEffects = function(trait){
  single_eff = read_csv(paste0("./data/", trait," traits/effectsSingleLocus.csv"))
  i5cM_eff   = read_csv(paste0("./data/", trait," traits/", trait,"_effectsInterval_5cM_mcmc.csv"))
  i10cM_eff  = read_csv(paste0("./data/", trait," traits/", trait,"_effectsInterval_10cM_mcmc.csv"))
  i15cM_eff  = read_csv(paste0("./data/", trait," traits/", trait,"_effectsInterval_15cM_mcmc.csv"))
  i20cM_eff  = read_csv(paste0("./data/", trait," traits/", trait,"_effectsInterval_20cM_mcmc.csv"))
  singleDIC  = read_csv(paste0("./data/", trait," traits/singleLocusDIC.csv"))
  intDIC     = read_csv(paste0("./data/", trait," traits/intervalMappingDIC.csv"))
 # nullDIC    = readRDS(file = paste0("./data/", trait," traits/nullDIC.rds"))
  single_eff = left_join(single_eff,
                         select(singleDIC, chrom, locus, singleLocus_DICDiff) %>%
                           rename(marker = locus),
                       by = c("chrom", "marker")) %>%
    rename(DICDiff  = singleLocus_DICDiff,
           ad_mean  = ad_post.mean, ad_lower = ad_l.95..CI, ad_upper = ad_u.95..CI,
           dm_mean  = dm_post.mean, dm_lower = dm_l.95..CI, dm_upper = dm_u.95..CI)
  i5cM_eff = left_join(i5cM_eff, select(intDIC, chrom, marker, DICDiff_5cM),
           by = c("chrom", "marker")) %>%
   rename(DICDiff  = DICDiff_5cM,
          ad_mean  = ad_post.mean, ad_lower = ad_l.95..CI, ad_upper = ad_u.95..CI,
          dm_mean  = dm_post.mean, dm_lower = dm_l.95..CI, dm_upper = dm_u.95..CI)
  i10cM_eff = left_join(i10cM_eff, select(intDIC, chrom, marker, DICDiff_10cM),
                       by = c("chrom", "marker")) %>%
    rename(DICDiff  = DICDiff_10cM,
           ad_mean  = ad_post.mean, ad_lower = ad_l.95..CI, ad_upper = ad_u.95..CI,
           dm_mean  = dm_post.mean, dm_lower = dm_l.95..CI, dm_upper = dm_u.95..CI)
  i15cM_eff = left_join(i15cM_eff, select(intDIC, chrom, marker, DICDiff_15cM),
                       by = c("chrom", "marker")) %>%
    rename(DICDiff  = DICDiff_15cM,
           ad_mean  = ad_post.mean, ad_lower = ad_l.95..CI, ad_upper = ad_u.95..CI,
           dm_mean  = dm_post.mean, dm_lower = dm_l.95..CI, dm_upper = dm_u.95..CI)
  i20cM_eff = left_join(i20cM_eff, select(intDIC, chrom, marker, DICDiff_20cM),
                       by = c("chrom", "marker")) %>%
    rename(DICDiff  = DICDiff_20cM,
           ad_mean  = ad_post.mean, ad_lower = ad_l.95..CI, ad_upper = ad_u.95..CI,
           dm_mean  = dm_post.mean, dm_lower = dm_l.95..CI, dm_upper = dm_u.95..CI)
  list(single = single_eff, i5cM = i5cM_eff, i10cM = i10cM_eff, i15cM = i15cM_eff,
       i20cM = i20cM_eff, singleDIC = singleDIC ,intDIC = intDIC
       #,nullDIC = nullDIC
       )
 }
trait_sets = c("area", "growth", "necropsy", "weight")
effects_list = llply(trait_sets, getEffects)
names(effects_list) = trait_sets

chrom_limits = 
	inner_join(select(effects_list[['necropsy']]$single, chrom, marker, count),
			   effects_list[[1]]$singleDIC %>% group_by(chrom) %>% summarise(marker = max(locus)),
           by = c("chrom", "marker")) %>%
	mutate(count = as.numeric(count) + 0.5) %>%
	 unique
chrom_limits$position = chrom_limits$count - c(31.5, diff(chrom_limits$count))/2

current_chrom = 6
x = effects_list[[4]]$i10cM
current_trait = "growth"
plotMainEffects = function(x, current_chrom, title = NULL, alpha = TRUE){
  plot_data =
    x %>%
    select(count, chrom, marker, trait,
           ad_mean, ad_lower, ad_upper,
           dm_mean, dm_lower, dm_upper, DICDiff) %>%
    filter(chrom == current_chrom) %>%
    rename(additive = ad_mean, dominance = dm_mean) %>%
    gather(type, value, additive, dominance) %>%
    mutate(trait = as.factor(eval(parse(text = paste0(current_trait, "_traits")))[trait]),
           chrom = as.factor(chrom))

    plot_data %>%
        filter(chrom == current_chrom) %>%
        ggplot(aes(marker, value)) +
        geom_point(size = 2) +
        geom_errorbar(data = filter(plot_data, type == "additive"),
                      aes(ymax = ad_upper, ymin = ad_lower), width = 0) +
        geom_errorbar(data = filter(plot_data, type == "dominance"),
                      aes(ymax = dm_upper, ymin = dm_lower), width = 0) +
        geom_hline(yintercept = 0) +
        facet_grid(trait ~ type, scales = "free") +
        ggtitle(paste(title))
}
plotMainEffects(effects_list[[2]]$single, 3, alpha = FALSE)

plotAllEffects = function(current_trait, prefix = "", ...){
    effect_types = c("i5cM","i10cM","i15cM","i20cM", "single")
    llply(effect_types, function(x) {
              for(chrom in 1:19){
                  save_plot(paste0("data/figures/", prefix, current_trait, "_",x,"_chrom", chrom, ".png"),
                            plotMainEffects(effects_list[[current_trait]][[x]],
                                            chrom, paste(current_trait, x, "chrom", chrom), ...),
                            base_aspect_ratio = 2, base_height = 8)
              }
                      })
}

# for(current_trait in trait_sets) plotAllEffects(current_trait, alpha = FALSE)
# for(current_trait in trait_sets) plotAllEffects(current_trait, "alpha_", alpha = TRUE)

interval_map_plots = list(
necropsy = effects_list[['necropsy']]$intDIC  %>%
	select(chrom, marker, contains("DICDiff")) %>%
	mutate(count = seq(353)) %>%
	rename( i5cM = DICDiff_5cM,
		   i10cM = DICDiff_10cM,
		   i15cM = DICDiff_15cM,
		   i20cM = DICDiff_20cM) %>%
	gather(interval, value, i20cM:i5cM) %>%
	mutate(interval = factor(interval, levels = c("i5cM", "i10cM", "i15cM", "i20cM")),
		   chrom = factor(chrom)) %>%
	filter(interval == "i10cM") %>%
	ggplot(aes(count, value)) +
	geom_line() + geom_hline(yintercept  = 15, linetype = "dashed") + geom_hline(yintercept = 0) +
	geom_vline(data = chrom_limits, linetype = "dotted", aes(xintercept = count)) +
	geom_text(data = chrom_limits, aes(x = position, y = 40, label = chrom)) + 
	labs(x = "marker", y = "DIC difference") + ggtitle("Necropsy organ weights") + annotate("text", x = 0, y  = 50, label = "Chromossome")
,
growth = effects_list[['growth']]$intDIC  %>%
	select(chrom, marker, contains("DICDiff")) %>%
	mutate(count = seq(353)) %>%
	rename( i5cM = DICDiff_5cM,
		   i10cM = DICDiff_10cM,
		   i15cM = DICDiff_15cM,
		   i20cM = DICDiff_20cM) %>%
	gather(interval, value, i20cM:i5cM) %>%
	mutate(interval = factor(interval, levels = c("i5cM", "i10cM", "i15cM", "i20cM")),
		   chrom = factor(chrom)) %>%
	filter(interval == "i10cM") %>%
	ggplot(aes(count, value)) +
	geom_line() + geom_hline(yintercept  = 15, linetype = "dashed") + geom_hline(yintercept = 0) +
	geom_vline(data = chrom_limits, linetype = "dotted", aes(xintercept = count)) +
	labs(x = "marker", y = "DIC difference") + ggtitle("Weekly growth")
,
area = effects_list[['area']]$intDIC  %>%
	select(chrom, marker, contains("DICDiff")) %>%
	mutate(count = seq(353)) %>%
	rename( i5cM = DICDiff_5cM,
		   i10cM = DICDiff_10cM,
		   i15cM = DICDiff_15cM,
		   i20cM = DICDiff_20cM) %>%
	gather(interval, value, i20cM:i5cM) %>%
	mutate(interval = factor(interval, levels = c("i5cM", "i10cM", "i15cM", "i20cM")),
		   chrom = factor(chrom)) %>%
	filter(interval == "i10cM") %>%
	ggplot(aes(count, value)) +
	geom_line() + geom_hline(yintercept  = 15, linetype = "dashed") + geom_hline(yintercept = 0) +
	geom_vline(data = chrom_limits, linetype = "dotted", aes(xintercept = count)) +
	labs(x = "marker", y = "DIC difference") + ggtitle("Mandible areas")
)
interval_DIC = plot_grid(interval_map_plots[[1]], 
		  interval_map_plots[[2]], 
		  interval_map_plots[[3]], ncol = 1)
save_plot("~/Dropbox/labbio/relatorios/fapesp/fapesp-relatorio-2016-10-30-BEPE/images/interval_DIC.png", interval_DIC, nrow = 3, base_height = 2.5, base_aspect_ratio = 5)

single_map_plots = list(
necropsy = effects_list[['necropsy']]$singleDIC  %>%
  select(chrom, locus, contains("DICDiff")) %>%
  mutate(count = seq(353)) %>%
  rename( single = singleLocus_DICDiff) %>%
  mutate(chrom = factor(chrom)) %>%
  ggplot(aes(count, single)) +
  geom_line() + geom_hline(yintercept  = 15, linetype = "dashed") + geom_hline(yintercept = 0) +
  geom_vline(data = chrom_limits, linetype = "dotted", aes(xintercept = count)) +
  geom_text(data = chrom_limits, aes(x = position, y = 110, label = chrom)) + 
  labs(x = "marker", y = "DIC difference") + ggtitle("Necropsy organ weights") + annotate("text", x = 0, y  = 125, label = "Chromossome"),
growth = effects_list[['growth']]$singleDIC  %>%
  select(chrom, locus, contains("DICDiff")) %>%
  mutate(count = seq(353)) %>%
  rename( single = singleLocus_DICDiff) %>%
  mutate(chrom = factor(chrom)) %>%
  ggplot(aes(count, single)) +
  geom_line() + geom_hline(yintercept  = 15, linetype = "dashed") + geom_hline(yintercept = 0) +
  geom_vline(data = chrom_limits, linetype = "dotted", aes(xintercept = count)) +
  geom_text(data = chrom_limits, aes(x = position, y = 110, label = chrom)) + 
  labs(x = "marker", y = "DIC difference") + ggtitle("Weekly growth"),
area = effects_list[['area']]$singleDIC  %>%
  select(chrom, locus, contains("DICDiff")) %>%
  mutate(count = seq(353)) %>%
  rename( single = singleLocus_DICDiff) %>%
  mutate(chrom = factor(chrom)) %>%
  ggplot(aes(count, single)) +
  geom_line() + geom_hline(yintercept  = 15, linetype = "dashed") + geom_hline(yintercept = 0) +
  geom_vline(data = chrom_limits, linetype = "dotted", aes(xintercept = count)) +
  geom_text(data = chrom_limits, aes(x = position, y = 110, label = chrom)) + 
  labs(x = "marker", y = "DIC difference") + ggtitle("Mandible areas")
)

single_DIC = plot_grid(single_map_plots[[1]], 
                         single_map_plots[[2]], 
                         single_map_plots[[3]], ncol = 1)
save_plot("~/Dropbox/labbio/relatorios/fapesp/fapesp-relatorio-2016-10-30-BEPE/images/regression_DIC.png", single_DIC, nrow = 3, base_height = 2.5, base_aspect_ratio = 5)

map_pos = read_csv("./data/markers/marker_positions.csv")[,2:3]
names(map_pos) <- c("SNP", "Pos")
detected_snps = 
  list(
    necropsy = cbind(map_pos, effects_list$necropsy$intDIC) %>% filter(DICDiff_10cM > 15) %>% select(chrom, marker, DICDiff_10cM, SNP, Pos),
    growth = cbind(map_pos, effects_list$growth$intDIC) %>% filter(DICDiff_10cM > 15) %>% select(chrom, marker, DICDiff_10cM, SNP, Pos),
area = cbind(map_pos, effects_list$area$intDIC) %>% filter(DICDiff_10cM > 15) %>% select(chrom, marker, DICDiff_10cM, SNP, Pos))

map_pos = read_csv("./data/markers/marker_positions.csv")[,2:3]

detected_snps$growth[-c(6, 9),]
predictGrowth <- function(x = detected_snps$growth, filename) {
  growth_effects = inner_join(effects_list$growth$single, x, by = c("chrom", "marker")) %>% 
    select(count, chrom, marker, trait, ad_mean) %>% spread(trait, ad_mean)
  
  LG = c(3.785,4.435,8.43,7.395,2.995,1.85,2.085)
  SM = c(3.31 ,2.98,3.82,2.175,0.765,1.165,0.51)
  F3 = sapply(growth_phen[,growth_traits], mean)
  
  SM_e = rowSums(-1 * t(as.matrix(growth_effects[,4:10]))) * sapply(growth_phen[,growth_traits], sd) + sapply(growth_phen[,growth_traits], mean)
  LG_e = rowSums(t(as.matrix(growth_effects[,4:10]))) * sapply(growth_phen[,growth_traits], sd) + sapply(growth_phen[,growth_traits], mean)
  
  growth_m = as.numeric(ddply(growth_phen, .(SEX), numcolwise(mean))[2,growth_traits])
  growth_f = as.numeric(ddply(growth_phen, .(SEX), numcolwise(mean))[1,growth_traits])
  
  growth_prediction = reshape2::melt(data.frame(trait = as.factor(growth_traits), 
                                                SM_Predicted = SM_e,
                                                SM_Observed = SM,
                                                LG_Predicted = LG_e,
                                                F3_Observed = F3,
                                                LG_Observed = LG)) %>% separate(variable, c("Line", "Type"))
  growth_pred_plot = ggplot(growth_prediction, aes(trait, value, group = interaction(Type, Line), color = Line, linetype = Type)) + geom_line(size = 1) + scale_x_discrete(labels = paste("Week", 1:7)) + labs(y = "Weekly growth (g)", x = "Start week")
  save_plot(paste0("~/Dropbox/labbio/relatorios/fapesp/fapesp-relatorio-2016-10-30-BEPE/images/", filename), growth_pred_plot, base_height = 4, base_aspect_ratio = 1.8)
  growth_pred_plot
}
predictGrowth(detected_snps$growth[-c(6, 9),], "growth_prediction.png")
predictGrowth(cbind(map_pos, effects_list$growth$intDIC) %>% filter(DICDiff_10cM > 10) %>% select(chrom, marker, DICDiff_10cM, SNP, Pos), "growth_prediction_2.png")
