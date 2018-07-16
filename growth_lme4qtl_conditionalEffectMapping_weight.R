setwd("/home/diogro/projects/mouse-qtls")
source('read_mouse_data.R')

#Rdatas_folder = "~/gdrive/LGSM_project_Rdatas/"
Rdatas_folder = "./Rdatas/"

ncores = 2
nthreads = 2
registerDoMC(ncores)
options(mc.cores = ncores)
setMKLthreads(nthreads)

growth_data = inner_join(growth_phen_std,
                        growth_markers,
                        by = "ID") 
growth_data = mutate(growth_data, FAMILY = as.factor(FAMILY))

markerMatrix = ldply(1:19, function(x) data.frame(chrom = x, marker = 1:loci_per_chrom[[x]]))
makeMarkerList = function(pos) paste(paste('chrom', pos[1],"_", c('A', 'D'), pos[2], sep = ''), collapse = ' + ')
markerList = alply(markerMatrix, 1, makeMarkerList)
markerPositions = cbind(markerMatrix, read_csv("./data/markers/marker_positions.csv")[,3])
names(markerPositions)[3] = "cM"

current_trait = 7
current_marker = 110
runInteractionModel <- function(current_trait, markers = 1:353){
  trait_term = weight_traits[1:current_trait]
  null_formula = paste(growth_traits[current_trait], " ~ (1|FAMILY)")
  makeInteractionMarkerList = function(pos){
    marker_term = paste('chrom', pos[1],"_", c('A', 'D'), pos[2], sep = '')
    paste(apply(expand.grid(trait_term, marker_term), 1, paste, collapse = ":"), collapse = " + ")
  }
  interactionMarkerList = alply(markerMatrix, 1, makeInteractionMarkerList)
  
  runInteractionModelSingleMaker = function(current_marker){
    direct_formula = paste(null_formula, markerList[current_marker], paste(trait_term, collapse = " + "), sep = " + ")
    direct_lmm = lmer(as.formula(direct_formula), data = growth_data, REML = FALSE)
    conditional_formula = paste(direct_formula, interactionMarkerList[current_marker], sep = " + ")
    conditional_lmm = lmer(as.formula(conditional_formula), 
                           data = growth_data, REML = FALSE)
    test = anova(direct_lmm, conditional_lmm)
    model_summary = summary(conditional_lmm)
    
    return(list(direct_lmm = direct_lmm,
                conditional_lmm = conditional_lmm, 
                test = test,
                p.value = test$'Pr(>Chisq)'[2],
                model_summary = model_summary))
  }
  llply(markers, runInteractionModelSingleMaker, .progress = "text")
}
conditional_gwas = vector("list", num_growth_traits)
names(conditional_gwas) = growth_traits
for(i in 1:(num_growth_traits)){
  print(growth_traits[i])
  conditional_gwas[[i]] = runInteractionModel(i)
}
save(conditional_gwas, file = "./Rdatas/growth_conditional_previous_weight.Rdata")
load("./Rdatas/growth_conditional_previous_weight.Rdata")

thresholds = read.csv("./data/F3_BONFERRONI_thresholds.csv")

Pvalues = function(i, ...){
    p.values = ldply(conditional_gwas[[i]], function(x) x$p.value) 
    p.values = rename(p.values, p_lrt = V1)
    p.values$significant = sapply(1:353, function(i) p.values$p_lrt[i] < thresholds[thresholds$chrom == markerMatrix[i, "chrom"], "p_value"])
    p.values$pos = markerPositions$cM
    p.values$snp = 1:353
    p.values$trait = growth_traits[i]
    p.values = cbind(markerMatrix, p.values)
    return(p.values)
}
p_values = ldply(1:7, Pvalues)
chrtable <- data.frame(table(filter(p_values, trait == "growth12")$chrom))
chrtable$Var1 <- as.character(chrtable$Var1)
chrtable <- chrtable[gtools:::mixedorder(chrtable$Var1), ]
oddchrom <- as.character(chrtable$Var1[seq(1, nrow(chrtable), 2)])
p_values$chrom_alt <- replace(p_values$chrom, p_values$chrom %in% oddchrom, 0)
p_values$chrom_alt <- replace(p_values$chrom_alt, p_values$chrom_alt != 0, 1)
dfmsplit <- split(filter(p_values, trait == "growth23"), filter(p_values, trait == "growth23")$chrom)
xbreaks <- sapply(dfmsplit, function(x) {
  midpoint <- length(x$snp)/2
  if (midpoint < 1) 
    midpoint <- 1
  return(x$snp[midpoint])
})

significantMarkerMatrix = read_csv("./data/growth_significant_markers.csv")
p_values_sig = inner_join(p_values, significantMarkerMatrix, by = c("chrom", "marker"))

(p1 = ggplot(p_values, aes(snp, -log10(p_lrt), color = as.factor(chrom_alt))) + 
  geom_point(aes(alpha = significant), size = 2) + 
  facet_wrap(~trait, ncol  = 1, scales = "free") + 
  scale_x_continuous(breaks = xbreaks) + labs(x = "Chromossome", y = "-log(p value)") + 
  guides(colour = FALSE) + geom_vline(data = p_values_sig, aes(xintercept = snp), color = "lightgrey") + 
  geom_hline(yintercept = -log10(thresholds[1,3]), size = 0.8, color = "grey", linetype = "dashed") + 
  theme(legend.position = "none") + scale_alpha_discrete(range = c(0.25, 1)) + scale_color_manual(values = c("black", "tomato3")))

save_plot("./data/growth_conditionalEffects_manhattan.png", p1, base_height = 7, base_aspect_ratio = 1.5)

significant_mask = vector("list", 6)
for(i in 1:6){
  aux = filter(p_values, significant == TRUE, trait == growth_traits[[i+1]])
  significant_mask[[i]] = cbind(aux$snp, aux$chrom, aux$marker, aux$p_lrt)
}
significant_mask[[6]]
x = list("1" = c(49, 67, 86, 202, 285, 296, 307),
         "2" = c(110, 156, 254, 277, 313, 344), 
         "3" = c(24, 56, 82, 105, 178, 210, 218, 267, 309),
         "4" = c(32, 51, 83, 89, 134, 196, 206, 220, 248, 320, 333),
         "5" = c(25, 116, 128, 156, 175, 221, 264, 302),
         "6" = c(67, 115, 125, 149, 167, 190, 228, 241, 283, 339))
for(i in 1:6){
  aux = significant_mask[[i]]
  significant_mask[[i]] = aux[aux[,1] %in% x[[i]],]
}
names(significant_mask) = growth_traits[-1]

conditional_effects = vector("list", 6)
for(trait in 1:6){
  conditional_effects[[trait]] = list()
  for(marker in 1:(nrow(significant_mask[[trait]]))){
      x = conditional_gwas[[trait]][significant_mask[[trait]][marker,1]][[1]]
      coef = x$model_summary$coefficients
      n_coef = nrow(coef)
      n_conditional = 2*trait
      get_lines = (n_coef - n_conditional + 1):n_coef
      coef = data.frame(coef[get_lines, c(1, 2, 5)])
      coef$trait = growth_traits[trait+1]
      coef$snp = significant_mask[[trait]][marker,1]
      coef$chrom = significant_mask[[trait]][marker,2]
      coef$marker = significant_mask[[trait]][marker,3]
      coef$id = paste(coef$chrom, coef$marker, sep = "_")
      coef$type = rep(c("additive", "dominance"), each = trait)
      coef$p_trait = rep(growth_traits[1:(trait)], 2)
      conditional_effects[[trait]][[marker]] = coef
      names(conditional_effects[[trait]][[marker]])[1:3] = c("mean", "sd", "p") 
  }
  conditional_effects[[trait]] = do.call(rbind, conditional_effects[[trait]])
}

conditional_effect_plots = vector("list", 6)
for (i in 1:6){
conditional_effect_plots[[i]]= ggplot(conditional_effects[[i]], aes(p_trait, mean, group = id, color = id)) + 
  geom_point(position = position_dodge(width = 0.3), size = 3) + 
  geom_pointrange(aes(ymin = mean - sd, ymax = mean + sd), position = position_dodge(width = 0.3)) +
  geom_hline(yintercept = 0) + 
  facet_wrap(~type) +
  scale_color_viridis_d(option ="B") +
  labs(x = "Week", y = "QTL effect") + scale_x_discrete(labels = 1:7)
}
all_conditiona = plot_grid(plotlist = conditional_effect_plots, labels = growth_traits[-1])
save_plot(here("data", "growth_conditional_effects.png"), all_conditiona, base_height = 5, base_aspect_ratio = 1.2, ncol = 3, nrow = 2)
