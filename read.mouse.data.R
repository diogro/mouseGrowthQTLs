library(gdata)
library(plyr)
library(dplyr)

raw.mouse_phen = read.csv("./data/F3phenotypes_uncorrected.csv", as.is = T)
raw.mouse_phen = select(raw.mouse_phen, c(ID, FATPAD:PAIR))

raw.mouse_phen = mutate(raw.mouse_phen,
                          grow12 = WEEK2 - WEEK1,
                          grow23 = WEEK3 - WEEK2,
                          grow34 = WEEK4 - WEEK3,
                          grow45 = WEEK5 - WEEK4,
                          grow56 = WEEK6 - WEEK5,
                          grow67 = WEEK7 - WEEK6,
                          grow78 = WEEK8 - WEEK7)

raw.mouse_meta = read.csv("./data/F3Phenotypes_further corrected family data_corrected litter sizes.csv", as.is = T)
names(raw.mouse_meta) = gsub('SexAN', 'SEX', names(raw.mouse_meta))
raw.mouse_meta = select(raw.mouse_meta, ID:COHORT)

raw.mouse_meta = raw.mouse_meta[raw.mouse_meta$ID %in% raw.mouse_phen$ID,]
raw.mouse_phen = raw.mouse_phen[raw.mouse_phen$ID %in% raw.mouse_meta$ID,]

mouse_phen = data.frame(arrange(raw.mouse_meta, ID), arrange(raw.mouse_phen, ID)[,-1])
raw.mouse_gen = llply(paste0("./data/genotypes/chrom", 1:19, ".csv"), read.csv, as.is = TRUE)
names(raw.mouse_gen) = paste0("chrom", 1:19)
mouse_gen = llply(raw.mouse_gen, function(x) x[x$ID %in% mouse_phen$ID,])
mouse_phen = mouse_phen[mouse_phen$ID %in% mouse_gen[[1]]$ID,]

rm(list = ls(pattern='raw'))
