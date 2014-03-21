library(plyr)
library(reshape2)
library(lme4)
library(ggplot2)

source('read.mouse.data.R')

m.data = melt(mouse.data, id.vars = names(mouse.data)[c(1:27, 35:127)])

ggplot(m.data, aes(variable, value, color = SEX)) + geom_boxplot()

null.formula = "value ~ 1 + variable * SEX + (1|FAMILY) + (0 + variable | FAMILY)"
mouse.model.no.gen = lmer(as.formula(null.formula), 
                          data = m.data, 
                          REML = FALSE)
G = VarCorr(mouse.model.no.gen)[[2]]

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
#save(all.loci, file= 'mouse.cromossome1.Rdata')
load("./mouse.cromossome1.Rdata")
significant = laply(all.loci, function(x) x$p.value < 0.05/31)

