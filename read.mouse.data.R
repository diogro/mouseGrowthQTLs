if(!require(install.load)) {install.packages("install.load"); library(install.load)}
install_load("plyr", "dplyr", "reshape2", "readr", "lme4")

##Growth traits

raw.growth_phen = read_csv("data/growth traits/F3Phenotypes_further corrected family data_corrected litter sizes.csv")
raw.growth_phen = tbl_df(select(raw.growth_phen, c(ID, FATPAD:LIVER, WEEK1:WEEK10)))

growth_phen = mutate(raw.growth_phen,
                          grow12 = WEEK2 - WEEK1,
                          grow23 = WEEK3 - WEEK2,
                          grow34 = WEEK4 - WEEK3,
                          grow45 = WEEK5 - WEEK4,
                          grow56 = WEEK6 - WEEK5,
                          grow67 = WEEK7 - WEEK6,
                          grow78 = WEEK8 - WEEK7)

raw.mouse_meta = read_csv("./data/growth traits/F3Phenotypes_further corrected family data_corrected litter sizes.csv")
names(raw.mouse_meta) = gsub('SexAN', 'SEX', names(raw.mouse_meta))
mouse_meta = select(raw.mouse_meta, ID:COHORT)
growth_phen = inner_join(mouse_meta, growth_phen, by = "ID")

markers = llply(paste0("./data/markers/chrom", 1:19, ".csv"), read_csv)
names(markers) = paste0("chrom", 1:19)

growth_phen = semi_join(growth_phen, markers[[1]], by = "ID")
growth_phen = arrange(growth_phen, ID)

growth_phen = select(growth_phen, ID, FAMILY, SEX, LSB, LSW, COHORT, grow12:grow78)
growth_phen = growth_phen[complete.cases(growth_phen),]
growth_markers = llply(markers, function(x) semi_join(x, growth_phen, by = "ID"))

traits = c( "grow12", "grow23", "grow34", "grow45", "grow56", "grow67", "grow78")
num_traits = length(traits)

m.data = melt(growth_phen, id.vars = names(growth_phen)[1:6])

null.formula = "value ~ 1 + variable * SEX + variable * LSB + variable * LSW + variable * COHORT"
mouse_no_fixed = lm(as.formula(null.formula), data = m.data)
m.data.std = m.data
m.data.std$value = residuals(mouse_no_fixed)

exclude = c(dim(m.data.std)[2]-1, dim(m.data.std)[2])
cast_formula = paste(paste(names(m.data.std[,-exclude]), collapse = " + "),
                     'variable',
                     sep = " ~ ")
growth_phen_std = tbl_df(dcast(m.data.std, as.formula(cast_formula)))

null.formula = "value ~ 1 + variable * SEX + variable * LSB + variable * LSW + variable * COHORT + (0 + variable|FAMILY)"
mouse_no_fixed = lmer(as.formula(null.formula), data = m.data)
m.data.std_mixed = m.data
m.data.std_mixed$value = residuals(mouse_no_fixed)

exclude = c(dim(m.data.std_mixed)[2]-1, dim(m.data.std_mixed)[2])
cast_formula = paste(paste(names(m.data.std[,-exclude]), collapse = " + "),
                     'variable',
                     sep = " ~ ")
growth_phen_std_mixed = tbl_df(dcast(m.data.std_mixed, as.formula(cast_formula)))

#Area traits

raw.area_phen = read_csv("data/area traits/areasF2F3.csv")
area_phen = inner_join(mouse_meta, raw.area_phen, by = "ID") %>% 
  semi_join(markers[[1]], by = "ID") %>% 
  select(ID, FAMILY, SEX, LSB, LSW, COHORT, area1:area7) %>%
  na.omit %>%
  arrange(ID)

area_markers = llply(markers, function(x) semi_join(x, area_phen, by = "ID"))

rm(list = ls(pattern='raw'))

