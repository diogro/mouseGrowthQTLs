library(gdata)
library(dplyr)

raw.mouse.phen = read.csv("./F3phenotypes_uncorrected.csv", as.is = T)
raw.mouse.phen = select(raw.mouse.phen, c(ID, FATPAD:PAIR))
  
raw.mouse.phen = mutate(raw.mouse.phen,
                          grow12 = WEEK2 - WEEK1,
                          grow23 = WEEK3 - WEEK2,
                          grow34 = WEEK4 - WEEK3,
                          grow45 = WEEK5 - WEEK4,
                          grow56 = WEEK6 - WEEK5,
                          grow67 = WEEK7 - WEEK6,
                          grow78 = WEEK8 - WEEK7)

raw.mouse.meta = read.csv("./F3Phenotypes_further corrected family data_corrected litter sizes.csv", as.is = T)
names(raw.mouse.meta) = gsub('SexAN', 'SEX', names(raw.mouse.meta))
raw.mouse.meta = select(raw.mouse.meta, ID:COHORT)

raw.mouse.meta = raw.mouse.meta[raw.mouse.meta$ID %in% raw.mouse.phen$ID,]
raw.mouse.phen = raw.mouse.phen[raw.mouse.phen$ID %in% raw.mouse.meta$ID,]

mouse.phen = data.frame(arrange(raw.mouse.meta, ID), arrange(raw.mouse.phen, ID)[,-1])
raw.mouse.gen = read.csv("./chrom1.csv")
mouse.gen = raw.mouse.gen[raw.mouse.gen$ID %in% mouse.phen$ID,]
mouse.phen = mouse.phen[mouse.phen$ID %in% mouse.gen$ID,]
mouse.data = cbind(arrange(mouse.phen, ID), arrange(mouse.gen, ID)[,-1])

rm('mouse.gen', 'mouse.phen', list = ls(pattern='raw'))
