setwd("/home/diogro/projects/mouse-qtls")
source('read_mouse_data.R')

#Rdatas_folder = "~/gdrive/LGSM_project_Rdatas/"
Rdatas_folder = "./Rdatas/"

install_load("doMC", "lme4qtl", "qvalue")
registerDoMC(10)

names(necropsy_phen_std)
necropsy_data = inner_join(necropsy_phen_std,
                        necropsy_markers,
                        by = "ID") %>%
  gather(variable, value, FATPAD:SPLEEN)
necropsy_data = mutate(necropsy_data, FAMILY = as.factor(FAMILY), variable = as.factor(variable))
null_formula = "value ~ variable + (0 + variable|FAMILY)"

makeMarkerList = function(pos) paste(paste('variable:chrom', pos[1],"_", c('A', 'D'), pos[2], sep = ''), collapse = ' + ')
markerMatrix = ldply(1:19, function(x) data.frame(chrom = x, marker = 1:loci_per_chrom[[x]]))
markerList = alply(markerMatrix, 1, makeMarkerList)
markerPositions = cbind(markerMatrix, read_csv("./data/markers/marker_positions.csv")[,3])
names(markerPositions)[3] = "cM"

runIntervalModel <- function(marker_term, null_formula){
    pos = na.omit(as.numeric(unlist(strsplit(unlist(as.character(marker_term)), "[^0-9]+"))))[3:2]
    focal_cm = dplyr::filter(markerPositions, chrom == pos[1], marker == pos[2])$cM
    possible_flank = filter(markerPositions, chrom == pos[1], abs(cM - focal_cm) > flank_dist)
    flanking = NULL
    for(line in 1:nrow(possible_flank)){
        if(possible_flank[line,"cM"] > focal_cm){ flanking = rbind(flanking, possible_flank[line,1:2]); break }
    }
    if(line > 1 & line != nrow(possible_flank)) flanking = rbind(possible_flank[line-1, 1:2], flanking)
    if(line == nrow(possible_flank)) flanking = rbind(possible_flank[line, 1:2], flanking)

    flanking_formula = paste(null_formula, paste(alply(flanking, 1, makeMarkerList), collapse = " + "), sep = " + ")

    flanking_model = relmatLmer(as.formula(flanking_formula),
                                data = necropsy_data,
                                REML = FALSE)


    genotype_formula = paste(flanking_formula, marker_term, sep = ' + ')

    focal_model = relmatLmer(as.formula(genotype_formula),
                             data = necropsy_data,
                             REML = FALSE)
    test = anova(focal_model, flanking_model)

    return(list(flanking = flanking_model,
                focal    = focal_model, 
                test = test,
                G = VarCorr(focal_model)[[1]],
                p.value = test$'Pr(>Chisq)'[2],
                model_summary = summary(focal_model)))
}
# flank_dist = 5
# model_file = paste0(Rdatas_folder, "necropsy_intervalMapping_", flank_dist, "cM.Rdata")
# intervalMapping = llply(markerList, runIntervalModel, null_formula, .parallel = TRUE)
# save(intervalMapping, file = model_file)
# 
# flank_dist = 10
# model_file = paste0(Rdatas_folder, "necropsy_intervalMapping_", flank_dist, "cM.Rdata")
# intervalMapping = llply(markerList, runIntervalModel, null_formula, .parallel = TRUE)
# save(intervalMapping, file = model_file)
# load(model_file)
# 
# flank_dist = 15
# model_file = paste0(Rdatas_folder, "necropsy_intervalMapping_", flank_dist, "cM.Rdata")
# intervalMapping = llply(markerList, runIntervalModel, null_formula, .parallel = TRUE)
# save(intervalMapping, file = model_file)
# 
# flank_dist = 20
# model_file = paste0(Rdatas_folder, "necropsy_intervalMapping_", flank_dist, "cM.Rdata")
# intervalMapping = llply(markerList, runIntervalModel, null_formula, .parallel = TRUE)
# save(intervalMapping, file = model_file)
# load(model_file)

all_effectsInterval = ldply(intervalMapping,
                            function(x){
                                n_traits = length(necropsy_traits)
                                coef_matrix =coef(x$model_summary) 
                                n_coef = nrow(coef_matrix)
                                ad = coef_matrix[(n_coef-2*n_traits + 1):(n_coef - n_traits),1]
                                ad_se = coef_matrix[(n_coef-2*n_traits + 1):(n_coef - n_traits),2]
                                dm = coef_matrix[(n_coef-n_traits + 1):n_coef,1]
                                dm_se = coef_matrix[(n_coef-n_traits + 1):n_coef,2]
                                data_frame(ad, ad_se, dm, dm_se)
                            }, .id = NULL, .parallel = TRUE) %>% tbl_df %>% mutate(count = rep(seq(intervalMapping), each = num_necropsy_traits)) %>% select(count, everything())

ldply(intervalMapping, function(x) -log10(x$p.value)) %>%
    filter(chrom == 1) %>%
    ggplot(aes(x =seq_along(V1), V1)) + geom_point() + geom_line()

effect_file = paste0("./data/necropsy traits/necropsy_effectsInterval_", flank_dist, "cM.csv")
write_csv(all_effectsInterval, effect_file)

thresholds = read.csv("./data/F3_BONFERRONI_thresholds.csv")

Pvalues = function(flank_dist, ...){
    model_file = paste0(Rdatas_folder, "necropsy_intervalMapping_", flank_dist, "cM.Rdata")
    load(model_file)
    p.values = ldply(intervalMapping, function(x) x$p.value) 
    p.values = rename(p.values, p_lrt = V1)
    p.values$significant = sapply(1:353, function(i) p.values$p_lrt[i] < thresholds[thresholds$chrom == markerMatrix[i, "chrom"], "p_value"])
    p.values$flank_dist = flank_dist
    p.values$pos = markerPositions$cM
    p.values$snp = 1:353
    return(p.values)
}
p_values = ldply(c(5, 10, 15, 20), Pvalues, .parallel = TRUE)
p_values$flank_dist_chr = factor(paste0("Flanking markers at ", p_values$flank_dist, "cM"), 
                                 levels = paste0("Flanking markers at ", c(20, 15, 10, 5), "cM"))
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

x = list("1" = c(5, 11, 29),
         "2" = c(10, 19, 24), 
         "3" = c(1, 25),
         "4" = c(12, 21),
         "6" = c(5, 13, 21),
         "7" = c(9, 13, 17),
         "8" = 2, 
         "9" = 6,
         "10"= c(4, 15),
         "11"= 12,
         "12"= c(3, 9),
         "13"= 9,
         "15"=1,
         "16"= 11,
         "17"= 5,
         "18"= 10, 
         "19"= 2)

significantMarkerMatrix = ldply(x, function(x) data.frame(marker = x), .id = "chrom")
write_csv(significantMarkerMatrix, "./data/necropsy_significant_markers.csv")

significantMarkerMatrix = ldply(x, function(x) data.frame(marker = x), .id = "chrom") 
significantMarkerMatrix$chrom = as.integer(as.character(significantMarkerMatrix$chrom))
p_values_sig = inner_join(p_values, significantMarkerMatrix, by = c("chrom", "marker"))
write_csv(significantMarkerMatrix, "./data/necropsy_significant_markers.csv")
filter(p_values, significant == TRUE, chrom == 4, flank_dist == 10)
library(ggman)
(p1 = ggplot(p_values, aes(snp, -log10(p_lrt), color = as.factor(chrom_alt))) + 
  geom_point(aes(alpha = significant), size = 2) + 
  facet_wrap(~flank_dist_chr, ncol  = 1, scales = "free") + 
  scale_x_continuous(breaks = xbreaks) + labs(x = "Chromossome", y = "-log(p value)") + 
  guides(colour = FALSE) + geom_vline(data = p_values_sig, aes(xintercept = snp), color = "lightgrey") + 
  geom_hline(yintercept = -log10(thresholds[1,3]), size = 0.8, color = "grey", linetype = "dashed") + 
  theme(legend.position = "none") + scale_alpha_discrete(range = c(0.25, 1)) + scale_color_manual(values = c("black", "tomato3")))
save_plot("./data/TalkStuff/necropsy_manhattan.png", p1, base_height = 7, base_aspect_ratio = 1.5)

significantMarkerList = alply(significantMarkerMatrix, 1, makeMarkerList)
significant_marker_term = paste(significantMarkerList, collapse = " + ")

