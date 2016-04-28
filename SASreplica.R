if(!require(plyr)) {install.packages("plyr"); library(plyr)}
if(!require(dplyr)) {install.packages("dplyr"); library(dplyr)}
if(!require(reshape2)) {install.packages("reshape2"); library(reshape2)}
if(!require(lme4)) {install.packages("lme4"); library(lme4)}
if(!require(ggplot2)) {install.packages("ggplot2"); library(ggplot2)}
if(!require(cowplot)) {install.packages("cowplot"); library(cowplot)}

source('read.mouse.data.R')

mouse.data = inner_join(select(mouse_phen, ID, FAMILY, SEX, LSB, LSW, COHORT, grow12:grow78),
                        mouse_gen[[1]], by = "ID")

m.data = melt(mouse.data, id.vars = names(mouse.data)[c(1:6, 14:106)])

data_plot = ggplot(m.data, aes(variable, value)) +
    geom_line(aes(variable, value, group = ID, color= SEX)) +
    geom_boxplot(aes(variable, value, group = variable)) + facet_wrap(~SEX)
save_plot("./data/data_test_plot.png", data_plot, base_height = 5, base_aspect_ratio = 2)

null.formula = "value ~ 1 + variable * SEX + variable * as.factor(LSB) + variable * as.factor(LSW) + variable * as.factor(COHORT) + (0 + variable|FAMILY)"
mouse.model.no.gen = lmer(as.formula(null.formula),
                          data = m.data,
                          REML = FALSE)
G.lme4 = VarCorr(mouse.model.no.gen)[[1]]
attr(G.lme4,"correlation") = NULL

runSingleLocusModel <- function(locus, null.formula){
    genotype.formula = paste(null.formula,
                             paste(paste('variable*', c('A', 'D', 'I'),
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
