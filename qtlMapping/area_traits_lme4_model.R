source("./read_mouse_data.R")

install_load("doMC")
registerDoMC(60)

area_data = inner_join(area_phen_std,
                       Reduce(inner_join, area_markers),
                       by = "ID") %>%
  gather(variable, value, area1:area7)

null_formula = "value ~ variable + (0 + variable|FAMILY)"
area_null_model = lmer(as.formula(null_formula),
                          data = area_data,
                          REML = FALSE)
G_lme4 = VarCorr(area_null_model)[[1]]
attr(G_lme4,"correlation") = NULL
(null_model_summary = summary(area_null_model))

makeMarkerList = function(pos) paste(paste('variable:chrom', pos[1],"_", c('A', 'D'), pos[2], sep = ''), collapse = ' + ')
markerMatrix = ldply(1:19, function(x) data.frame(chrom = x, locus = 1:loci_per_chrom[[x]]))
markerList = alply(markerMatrix, 1, makeMarkerList)

#marker_term = markerList[[50]]
runSingleLocusModel <- function(marker_term, null_formula){
    genotype_formula = paste(null_formula, marker_term, sep = ' + ')
    single_locus_model = lmer(as.formula(genotype_formula),
                                    data = area_data,
                                    REML = FALSE)
    coef(lme4_summary )
    coef(lmerTest_summary)
    test = lmerTest::anova(area_null_model, single_locus_model)
    return(list(model = single_locus_model,
                model_summary = summary(single_locus_model),
                test = test,
                G = VarCorr(single_locus_model)[[1]],
                p.value = test$'Pr(>Chisq)'[2]))
}
all_loci = llply(markerList, runSingleLocusModel, null_formula, .parallel = TRUE)
#save(all_loci, file = "./data/Rdatas/all_loci_lme4.Rdata")
all_effects = ldply(all_loci, function(x) coef(x$model_summary)[8:21,1:2])
names(all_effects)
all_effects %>% filter(chrom == 2, locus == 9)

