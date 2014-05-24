library(plyr)
library(dplyr)
library(reshape2)
library(lme4)
library(ggplot2)

source('read.mouse.data.R')

mouse.data = select(mouse.data, ID, FAMILY, SEX, LSB, LSW, COHORT, grow12:grow78, A1:I31)
mouse.data = mouse.data[complete.cases(mouse.data),]

m.data = melt(mouse.data, id.vars = names(mouse.data)[c(1:6, 14:106)])

ggplot(m.data, aes(variable, value, color = SEX)) + geom_boxplot() + facet_wrap(~SEX)

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
    print(test)

    return(list(model = mouse.model, anova = test, G = VarCorr(mouse.model.no.gen)[[2]], p.value = test$'Pr(>Chisq)'[2]))
}

#all.loci = alply(1:31, 1, runSingleLocusModel, null.formula)
#save(all.loci, file= './Rdatas/mouse.cromossome1.Rdata')
load("./Rdatas/mouse.cromossome1.Rdata")
significant = laply(all.loci, function(x) x$p.value < 0.05/31)
