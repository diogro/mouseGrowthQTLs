setwd("/home/diogro/projects/mouse-qtls")
source('read_mouse_data.R')

#Rdatas_folder = "~/gdrive/LGSM_project_Rdatas/"
Rdatas_folder = "./Rdatas/"

ncores = 1
nthreads = 8
registerDoMC(ncores)
options(mc.cores = ncores)
setMKLthreads(nthreads)

growth_data = inner_join(growth_phen_std,
                        growth_markers,
                        by = "ID") 
growth_data = mutate(growth_data, FAMILY = as.factor(FAMILY), variable = as.factor(variable))

markerMatrix = ldply(1:19, function(x) data.frame(chrom = x, marker = 1:loci_per_chrom[[x]]))
makeMarkerList = function(pos) paste(paste('chrom', pos[1],"_", c('A', 'D'), pos[2], sep = ''), collapse = ' + ')
markerList = alply(markerMatrix, 1, makeMarkerList)
markerPositions = cbind(markerMatrix, read_csv("./data/markers/marker_positions.csv")[,3])
names(markerPositions)[3] = "cM"

current_trait = 2
runInteractionModel <- function(current_trait, markers = 1:353){
  trait_term = growth_traits[1:(current_trait-1)]
  null_formula = paste(growth_traits[current_trait], " ~ (1|FAMILY)")
  makeInteractionMarkerList = function(pos){
    marker_term = paste('chrom', pos[1],"_", c('A', 'D'), pos[2], sep = '')
    paste(apply(expand.grid(trait_term, marker_term), 1, paste, collapse = ":"), collapse = " + ")
  }
  interactionMarkerList = alply(markerMatrix, 1, makeInteractionMarkerList)
  
  runInteractionModelSingleMaker = function(current_marker){
    direct_formula = paste(null_formula, markerList[current_marker], trait_term, sep = " + ")
    direct_lmm = lmer(as.formula(direct_formula), data = growth_data, REML = FALSE)
    conditional_formula = paste(direct_formula, interactionMarkerList[current_marker], sep = " + ")
    conditional_lmm = lmer(growth23 ~ (1|FAMILY) + chrom1_A1 + chrom1_D1 + 
                             growth12:chrom1_A1 + growth12:chrom1_D1, 
                           data = growth_data, REML = FALSE)
    test = anova(direct_lmm, conditional_lmm)
    
    return(list(direct_lmm = direct_lmm,
                conditional_lmm = conditional_lmm, 
                test = test,
                p.value = test$'Pr(>Chisq)'[2],
                model_summary = summary(conditional_lmm)))
  }
  llply(markers, runInteractionModelSingleMaker, .progress = "text")
}
# conditional_gwas = vector("list", num_growth_traits-1)
# names(conditional_gwas) = growth_traits[-1]
# for(i in 1:(num_growth_traits-1)){
#   print(growth_traits[i+1])
#   conditional_gwas[[i]] = runInteractionModel(i+1)
# }
# save(conditional_gwas, file = "./Rdatas/growth_conditional_all_previous.Rdata")
load("./Rdatas/growth_conditional_all_previous.Rdata")

thresholds = read.csv("./data/F3_BONFERRONI_thresholds.csv")

Pvalues = function(i, ...){
    p.values = ldply(conditional_gwas[[i]], function(x) x$p.value) 
    p.values = rename(p.values, p_lrt = V1)
    p.values$significant = sapply(1:353, function(i) p.values$p_lrt[i] < thresholds[thresholds$chrom == markerMatrix[i, "chrom"], "p_value"])
    p.values$pos = markerPositions$cM
    p.values$snp = 1:353
    p.values$trait = growth_traits[i+1]
    p.values = cbind(markerMatrix, p.values)
    return(p.values)
}
p_values = ldply(1:6, Pvalues)
chrtable <- data.frame(table(p_values$chrom))
chrtable$Var1 <- as.character(chrtable$Var1)
chrtable <- chrtable[gtools:::mixedorder(chrtable$Var1), ]
oddchrom <- as.character(chrtable$Var1[seq(1, nrow(chrtable), 2)])
p_values$chrom_alt <- replace(p_values$chrom, p_values$chrom %in% oddchrom, 0)
p_values$chrom_alt <- replace(p_values$chrom_alt, p_values$chrom_alt != 0, 1)
dfmsplit <- split(p_values, p_values$chrom)
xbreaks <- sapply(dfmsplit, function(x) {
  midpoint <- length(x$snp)/8
  if (midpoint < 1) 
    midpoint <- 1
  return(x$snp[midpoint])
})

x = list("1" = c(4, 19, 29),
         "2" = c(23), 
         "4" = c(15, 21),
         "5" = c(11, 16, 18),
         "6" = c(4, 14, 18, 22),
         "7" = 9,
         "8" = c(2, 8, 12),
         "9" = 4,
         "10"= c(4, 11, 16),
         "11"= c(14),
         "12"= c(15, 19),
         "13"= 6,
         "14"= c(3, 8, 14),
         "15"= 15,
         "17"= 5,
         "18"= c(5, 12))
length(unlist(x))

significantMarkerMatrix = ldply(x, function(x) data.frame(marker = x), .id = "chrom") 
significantMarkerMatrix$chrom = as.integer(as.character(significantMarkerMatrix$chrom))
p_values_sig = inner_join(p_values, significantMarkerMatrix, by = c("chrom", "marker"))
#write_csv(significantMarkerMatrix, "./data/growth_significant_markers.csv")
library(ggman)
(p1 = ggplot(p_values, aes(snp, -log10(p_lrt), color = as.factor(chrom_alt))) + 
  geom_point(aes(alpha = significant), size = 2) + 
  facet_wrap(~trait, ncol  = 1, scales = "free") + 
  scale_x_continuous(breaks = xbreaks) + labs(x = "Chromossome", y = "-log(p value)") + 
  guides(colour = FALSE) + geom_vline(data = p_values_sig, aes(xintercept = snp), color = "lightgrey") + 
  geom_hline(yintercept = -log10(thresholds[1,3]), size = 0.8, color = "grey", linetype = "dashed") + 
  theme(legend.position = "none") + scale_alpha_discrete(range = c(0.25, 1)) + scale_color_manual(values = c("black", "tomato3")))

#save_plot("./data/TalkStuff/growth_manhattan.png", p1, base_height = 7, base_aspect_ratio = 1.5)

significantMarkerList = alply(significantMarkerMatrix, 1, makeMarkerList)
significant_marker_term = paste(significantMarkerList, collapse = " + ")

