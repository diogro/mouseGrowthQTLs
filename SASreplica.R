library(gdata)
library(dplyr)
library(reshape2)
library(lme4)

raw.mouse.phen = read.xls("./F3Phenotypes_further corrected family data_corrected litter sizes.xls")
names(raw.mouse.phen) = gsub('SexAN', 'SEX', names(raw.mouse.phen))

raw.mouse.phen = mutate(raw.mouse.phen,
                        grow12 = WEEK2 - WEEK1,
                        grow23 = WEEK3 - WEEK2,
                        grow34 = WEEK4 - WEEK3,
                        grow45 = WEEK5 - WEEK4,
                        grow56 = WEEK6 - WEEK5,
                        grow67 = WEEK7 - WEEK6,
                        grow78 = WEEK8 - WEEK7)

mouse.phen = select(raw.mouse.phen, ID:LIVER, grow12:grow78)

raw.mouse.gen = read.xls("./chrom1.xls")
mouse.gen = raw.mouse.gen[raw.mouse.gen$ID %in% mouse.phen$ID,]
mouse.phen = mouse.phen[mouse.phen$ID %in% mouse.gen$ID,]
mouse.data = cbind(arrange(mouse.phen, ID), arrange(mouse.gen, ID)[,-1])

m.data = melt(mouse.data, id.vars = names(mouse.data)[c(1:27, 41:133)])

null.formula = "value ~ variable * SEX + (1|FAMILY)"
locus = 37
genotype.formula = paste(null.formula, 
                         paste(paste(c('A', 'D', 'I'), 
                                     locus, sep = ''), collapse = ' + '), 
                         sep = ' + ')
mouse.model.no.gen = lmer(as.formula(null.formula), 
                          data = m.data, 
                          REML=FALSE)
mouse.model = lmer(as.formula(genotype.formula), 
                   data = m.data, 
                   REML=FALSE)
anova(mouse.model.no.gen, mouse.model)

