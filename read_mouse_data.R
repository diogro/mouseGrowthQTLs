if(!require(install.load)) {install.packages("install.load"); library(install.load)}
install_load("plyr", "dplyr", "tidyr", "readr", "lme4")

## Meta data

raw.mouse_meta = read_csv("./data/growth traits/F3Phenotypes_further corrected family data_corrected litter sizes.csv")
names(raw.mouse_meta) = gsub('SexAN', 'SEX', names(raw.mouse_meta))
mouse_meta = select(raw.mouse_meta, ID:COHORT)

## Markers

markers = llply(paste0("./data/markers/chrom", 1:19, ".csv"), read_csv)
names(markers) = paste0("chrom", 1:19)

##Growth traits

raw.growth_phen = read_csv("data/growth traits/F3Phenotypes_further corrected family data_corrected litter sizes.csv")
raw.growth_phen = tbl_df(select(raw.growth_phen, c(ID, WEEK1:WEEK10)))

growth_phen = inner_join(mouse_meta, raw.growth_phen, by = "ID") %>% 
  semi_join(markers[[1]], by = "ID") %>% 
  mutate(grow12 = WEEK2 - WEEK1,
         grow23 = WEEK3 - WEEK2,
         grow34 = WEEK4 - WEEK3,
         grow45 = WEEK5 - WEEK4,
         grow56 = WEEK6 - WEEK5,
         grow67 = WEEK7 - WEEK6,
         grow78 = WEEK8 - WEEK7) %>%
  select(ID, FAMILY, SEX, LSB, LSW, COHORT, grow12:grow78) %>%
  na.omit %>%
  arrange(ID)

growth_markers = llply(markers, function(x) semi_join(x, growth_phen, by = "ID"))

growth_traits = c("grow12", "grow23", "grow34", "grow45", "grow56", "grow67", "grow78")
num_growth_traits = length(growth_traits)

m_growth_phen = gather(growth_phen, variable, value, grow12:grow78)

null.formula = "value ~ 1 + variable * SEX + variable * LSB + variable * LSW + variable * COHORT"
mouse_no_fixed = lm(as.formula(null.formula), data = m_growth_phen)
m_growth_phen_std = m_growth_phen
m_growth_phen_std$value = residuals(mouse_no_fixed)

growth_phen_std = spread(m_growth_phen_std, variable, value)

#Area traits

raw.area_phen = read_csv("data/area traits/areasF2F3.csv")
area_phen = inner_join(mouse_meta, raw.area_phen, by = "ID") %>% 
  semi_join(markers[[1]], by = "ID") %>% 
  select(ID, FAMILY, SEX, LSB, LSW, COHORT, area1:area7) %>%
  na.omit %>%
  arrange(ID)

area_markers = llply(markers, function(x) semi_join(x, area_phen, by = "ID"))

area_traits = paste0("area", 1:7)
num_area_traits = length(area_traits)

m_area_phen = gather(area_phen, variable, value, area1:area7)

null.formula = "value ~ 1 + variable * SEX + variable * LSB + variable * LSW + variable * COHORT"
mouse_no_fixed = lm(as.formula(null.formula), data = m_area_phen)
m_area_phen_std = m_area_phen
m_area_phen_std$value = residuals(mouse_no_fixed)

area_phen_std = spread(m_area_phen_std, variable, value)

rm(list = ls(pattern='raw'))
rm(list = c("null.formula", 'mouse_no_fixed'))