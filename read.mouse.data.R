if(!require(gdata)) {install.packages("gdata"); library(gdata)}
if(!require(plyr)) {install.packages("plyr"); library(plyr)}
if(!require(dplyr)) {install.packages("dplyr"); library(dplyr)}
if(!require(reshape2)) {install.packages("reshape2"); library(reshape2)}
if(!require(readr)) {install.packages("readr"); library(readr)}

raw.mouse_phen = read_csv("./data/F3phenotypes_uncorrected.csv")
raw.mouse_phen = tbl_df(select(raw.mouse_phen, c(ID, FATPAD:LIVER, WEEK1:WEEK10)))

raw.mouse_phen = mutate(raw.mouse_phen,
                          grow12 = WEEK2 - WEEK1,
                          grow23 = WEEK3 - WEEK2,
                          grow34 = WEEK4 - WEEK3,
                          grow45 = WEEK5 - WEEK4,
                          grow56 = WEEK6 - WEEK5,
                          grow67 = WEEK7 - WEEK6,
                          grow78 = WEEK8 - WEEK7)

raw.mouse_meta = read_csv("./data/F3Phenotypes_further corrected family data_corrected litter sizes.csv")
names(raw.mouse_meta) = gsub('SexAN', 'SEX', names(raw.mouse_meta))
raw.mouse_meta = select(raw.mouse_meta, ID:COHORT)
mouse_phen = inner_join(raw.mouse_meta, raw.mouse_phen, by = "ID")

raw.mouse_gen = llply(paste0("./data/genotypes/chrom", 1:19, ".csv"), read_csv)
names(raw.mouse_gen) = paste0("chrom", 1:19)
mouse_gen = llply(raw.mouse_gen, function(x) semi_join(x, mouse_phen, by = "ID"))
mouse_gen = llply(mouse_gen, function(x) arrange(x, ID))

mouse_phen = semi_join(mouse_phen, mouse_gen[[1]], by = "ID")
mouse_phen = arrange(mouse_phen, ID)

rm(list = ls(pattern='raw'))

mouse_phen = select(mouse_phen, ID, FAMILY, SEX, LSB, LSW, COHORT, grow12:grow78)
complete_rows = complete.cases(mouse_phen)
mouse_phen = mouse_phen[complete_rows,]
mouse_gen = llply(mouse_gen, function(x) x[complete_rows,])

num_traits = 7
traits = c( "grow12", "grow23", "grow34", "grow45", "grow56", "grow67", "grow78")

m.data = melt(mouse_phen, id.vars = names(mouse_phen)[1:6])

null.formula = "value ~ 1 + variable * SEX + variable * LSB + variable * LSW + variable * COHORT"
mouse_no_fixed = lm(as.formula(null.formula), data = m.data)
m.data.std = m.data
m.data.std$value = residuals(mouse_no_fixed)

exclude = c(dim(m.data.std)[2]-1, dim(m.data.std)[2])
cast_formula = paste(paste(names(m.data.std[,-exclude]), collapse = " + "),
                     'variable',
                     sep = " ~ ")
mouse_phen_std = tbl_df(dcast(m.data.std, as.formula(cast_formula)))
