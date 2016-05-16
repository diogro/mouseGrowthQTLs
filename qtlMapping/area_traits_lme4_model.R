source("./read_mouse_data.R")

null.formula = "value ~ 1 + (0 + variable|FAMILY)"
mouse.model.no.gen = lmer(as.formula(null.formula),
                          data = m_area_phen_std,
                          REML = FALSE)
G_lme4 = VarCorr(mouse.model.no.gen)[[1]]
attr(G_lme4,"correlation") = NULL

runSingleLocusModel <- function(locus, null.formula){
    genotype.formula = paste(null.formula,
                             paste(paste('variable*', c('A', 'D'),
                                         locus, sep = ''), collapse = ' + '),
                             sep = ' + ')
    mouse.model = lmer(as.formula(genotype.formula),
                       data = m.data,
                       REML = FALSE)
    test = anova(mouse.model.no.gen, mouse.model)
    print(locus)
    return(list(model = mouse.model,
                anova = test,
                G = VarCorr(mouse.model.no.gen)[[1]],
                p.value = test$'Pr(>Chisq)'[2]))
}
names(m.data)
all.loci = alply(1:31, 1, runSingleLocusModel, null.formula, .parallel = TRUE)
save(all.loci, file= './data/Rdatas/mouse.cromossome1.Rdata')
#load("./data/Rdatas/mouse.cromossome1.Rdata")
significant = laply(all.loci, function(x) x$p.value < 0.05/31)
pvalue_plot = ldply(all.loci, function(x) c(x$p.value < 0.05/31, log(x$p.value))) %>%
    ggplot(aes(X1, V2, fill = V1)) + geom_bar(stat = "identity")
save_plot("./data/pvalue_test_plot.png", pvalue_plot, base_height = 5)
