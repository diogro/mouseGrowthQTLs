source("./read_mouse_data.R")

install_load("doMC")
registerDoMC(2)


area_data = inner_join(area_phen_std, 
                       Reduce(inner_join, area_markers), 
                       by = "ID") %>%
  gather(variable, value, area1:area7)

null_formula = "value ~ 1 + (0 + variable|FAMILY)"
area_null_model = lmer(as.formula(null_formula),
                          data = area_data,
                          REML = FALSE)
G_lme4 = VarCorr(area_null_model)[[1]]
attr(G_lme4,"correlation") = NULL
summary(area_null_model)

chrom = 10
locus = 1
runSingleLocusModel <- function(locus, chrom, null_formula){
    genotype_formula = paste(null_formula,
                             paste(paste('variable*chrom',chrom,"_", c('A', 'D'),
                                         locus, sep = ''), collapse = ' + '),
                             sep = ' + ')
    single_locus_model = lmer(as.formula(genotype_formula),
                       data = area_data,
                       REML = FALSE)
    test = lmerTest::anova(area_null_model, single_locus_model)
    return(list(model = single_locus_model,
                model_summary = summary(single_locus_model),
                test = test,
                G = VarCorr(single_locus_model)[[1]],
                p.value = test$'Pr(>Chisq)'[2]))
}
x = runSingleLocusModel(1, 10, null_formula)
x$model_summary
all_loci = alply(seq(chroms), 1, 
                 function(chrom) alply(seq(loci_per_chrom[chrom]), 1, 
                                       runSingleLocusModel, chrom, null_formula), 
                 .parallel = TRUE)
names(all_loci) = names(markers)
