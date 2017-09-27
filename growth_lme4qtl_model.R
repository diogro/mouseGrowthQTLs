setwd("/home/diogro/projects/mouse-qtls/")
source("./read_mouse_data.R")

install_load("doMC", "lme4qtl")
registerDoMC(10)

growth_data = inner_join(growth_phen_std,
                        growth_markers,
                        by = "ID") %>%
  gather(variable, value, growth12:growth78)
growth_data = mutate(growth_data, FAMILY = as.factor(FAMILY), variable = as.factor(variable))
null_formula = "value ~ variable + (0 + variable|FAMILY)"

growth_null_model = relmatLmer(as.formula(null_formula),
                               data = growth_data,
                               REML = FALSE)
summary(growth_null_model)
G_ML = VarCorr(growth_null_model)[[1]]
attr(G_ML,"correlation") = NULL
(null_model_summary = summary(growth_null_model))

makeMarkerList = function(pos) paste(paste('variable:chrom', pos[1],"_", c('A', 'D'), pos[2], sep = ''), collapse = ' + ')
markerMatrix = ldply(1:19, function(x) data.frame(chrom = x, locus = 1:loci_per_chrom[[x]]))
markerList = alply(markerMatrix, 1, makeMarkerList)

runSingleLocusModel <- function(marker_term, null_formula){
    genotype_formula = paste(null_formula, marker_term, sep = ' + ')
    single_locus_model = relmatLmer(as.formula(genotype_formula),
                                    data = growth_data,
                                    REML = FALSE)
    test = anova(single_locus_model, growth_null_model)
    return(list(model = single_locus_model,
                model_summary = summary(single_locus_model),
                test = test,
                G = VarCorr(single_locus_model)[[1]],
                p.value = test$'Pr(>Chisq)'[2]))
}
all_loci = llply(markerList, runSingleLocusModel, null_formula, .parallel = TRUE)
save(all_loci, file = "./Rdatas/all_loci_lme4qtl.Rdata")
load("./Rdatas/all_loci_lme4qtl.Rdata")
x = all_loci[[25]]
x$p.value 
all_effects = ldply(all_loci,
                    function(x){
                        ad = coef(x$model_summary)[8:14,1]
                        ad_se = coef(x$model_summary)[8:14,2]
                        #p_ad = coef(x$model_summary)[8:14,5]
                        dm = coef(x$model_summary)[15:21,1]
                        dm_se = coef(x$model_summary)[15:21,2]
                        #p_dm = coef(x$model_summary)[15:21,5]
                        data_frame(ad, ad_se, dm, dm_se)
                    }, .inform = TRUE ) %>% mutate(count = 1:nrow(.))
names(all_effects)
write_csv(all_effects,"growht_effects.csv" )
all_effects %>% filter(chrom == 2, locus == 9)
library(tidyverse)
ldply(all_loci, function(x) -log10(x$p.value)) %>% ggplot(aes(x =seq_along(V1), V1, color = as.factor(chrom))) + geom_point()


