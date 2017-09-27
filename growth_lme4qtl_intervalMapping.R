setwd("/home/diogro/projects/mouse-qtls")
source('read_mouse_data.R')

#Rdatas_folder = "~/gdrive/LGSM_project_Rdatas/"
Rdatas_folder = "./Rdatas/"

install_load("doMC", "lme4qtl", "qvalue")
registerDoMC(10)

growth_data = inner_join(growth_phen_std,
                        growth_markers,
                        by = "ID") %>%
  gather(variable, value, growth12:growth78)
growth_data = mutate(growth_data, FAMILY = as.factor(FAMILY), variable = as.factor(variable))
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
                                data = growth_data,
                                REML = FALSE)


    genotype_formula = paste(flanking_formula, marker_term, sep = ' + ')

    focal_model = relmatLmer(as.formula(genotype_formula),
                             data = growth_data,
                             REML = FALSE)
    test = anova(focal_model, flanking_model)

    return(list(flanking = flanking_model,
                focal    = focal_model, 
                test = test,
                G = VarCorr(focal_model)[[1]],
                p.value = test$'Pr(>Chisq)'[2],
                model_summary = summary(focal_model)))
}
flank_dist = 10
model_file = paste0(Rdatas_folder, "growth_intervalMapping_", flank_dist, "cM.Rdata")
#intervalMapping = llply(markerList, runIntervalModel, null_formula, .parallel = TRUE)
#save(intervalMapping, file = model_file)
load(model_file)

x = intervalMapping[[10]]
coef(x$model_summary)
all_effectsInterval = ldply(intervalMapping,
                            function(x){
                                n_traits = length(growth_traits)
                                coef_matrix =coef(x$model_summary) 
                                n_coef = nrow(coef_matrix)
                                ad = coef_matrix[(n_coef-2*n_traits + 1):(n_coef - n_traits),1]
                                ad_se = coef_matrix[(n_coef-2*n_traits + 1):(n_coef - n_traits),2]
                                dm = coef_matrix[(n_coef-n_traits + 1):n_coef,1]
                                dm_se = coef_matrix[(n_coef-n_traits + 1):n_coef,2]
                                data_frame(ad, ad_se, dm, dm_se)
                            }, .id = NULL, .parallel = TRUE) %>% tbl_df %>% mutate(count = rep(seq(intervalMapping), each = num_growth_traits)) %>% select(count, everything())

ldply(intervalMapping, function(x) -log10(x$p.value)) %>%
    filter(chrom == 1) %>%
    ggplot(aes(x =seq_along(V1), V1)) + geom_point() + geom_line()

effect_file = paste0("./data/growth traits/growth_effectsInterval_", flank_dist, "cM.csv")
write_csv(all_effectsInterval, effect_file)
dmSend(paste("Finished growth interval mapping with flanking", flank_dist, "cM"), "diogro")

hist(p.values$V1, nclass = 20)
summary(qobj)
plot(qobj)
qobj$significant


Pvalues = function(flank_dist, ...){
    model_file = paste0(Rdatas_folder, "growth_intervalMapping_", flank_dist, "cM.Rdata")
    load(model_file)
    p.values = ldply(intervalMapping, function(x) x$p.value) 
    qobj = qvalue(p.values$V1, ...)
    p.values$q_values = qobj$qvalues
    p.values$significant = qobj$significant
    p.values$flank_dist = flank_dist
    return(p.values)
}
p_values = ldply(c(5, 10, 15, 20), Pvalues, 0.05)
p_values %>% 
    filter(flank_dist == 5) %>% 
    ggplot(aes(x =seq_along(V1), V1, color = as.factor(chrom))) + geom_point() + geom_line()

p_values %>%
    select(chrom, marker, significant, flank_dist) %>%
    filter(chrom == 19) %>%
    spread(flank_dist, significant)

ldply(intervalMapping, function(x) -log10(x$p.value)) %>%
    filter(chrom == 19) %>%
    ggplot(aes(x =seq_along(V1), V1)) + geom_point() + geom_line()

x = list("1" = c(4, 19, 29),
         "2" = 23, 
         "4" = 15,
         "5" = 11,
         "6" = c(4, 14, 19),
         "7" = 9,
         "8" = c(2, 8, 12),
         "9" = 5,
         "10"= c(4, 11, 16),
         "11"= c(14, 17),
         "12"= c(15, 19),
         "14"= c(3, 8, 14),
         "15"= 15,
         "17"= 5,
         "18"= c(4, 12))

significantMarkerMatrix = ldply(x, function(x) data.frame(marker = x))
