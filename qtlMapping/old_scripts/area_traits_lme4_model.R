setwd("/home/diogro/projects/mouse-qtls/")
source("./read_mouse_data.R")

install_load("doMC", "VCA")
registerDoMC(4)

#area_data = inner_join(area_phen_std,
                       #area_markers,
                       #by = "ID") %>%
  #gather(variable, value, area1:area7)
area_data = area_phen_std %>%
  gather(variable, value, area1:area7)

null_formula = "value ~ variable + (0 + variable|FAMILY)"
area_null_model = lmer(as.formula(null_formula),
                          data = area_data,
                          REML = FALSE)
summary(area_null_model)
G_ML = VarCorr(area_null_model)[[1]]
attr(G_ML,"correlation") = NULL
(null_model_summary = summary(area_null_model))

area_data = mutate(area_data, FAMILY = as.factor(FAMILY), variable = as.factor(variable))
area_vcf = remlMM(value ~ variable + variable:(FAMILY), as.data.frame(area_data))


source("./qtlMapping/check_G.R")

makeMarkerList = function(pos) paste(paste('variable:chrom', pos[1],"_", c('A', 'D'), pos[2], sep = ''), collapse = ' + ')
markerMatrix = ldply(1:19, function(x) data.frame(chrom = x, locus = 1:loci_per_chrom[[x]]))
markerList = alply(markerMatrix, 1, makeMarkerList)

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
#all_loci = llply(markerList, runSingleLocusModel, null_formula, .parallel = TRUE)
#save(all_loci, file = "./data/Rdatas/all_loci_lme4.Rdata")
load("./data/Rdatas/all_loci_lme4.Rdata")
x = all_loci[[25]]
x$model_summary
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
all_effects %>% filter(chrom == 2, locus == 9)
