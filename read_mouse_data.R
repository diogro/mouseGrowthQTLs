if(!require(rstan)){install.packages("rstan"); library(rstan)}
if(!require(plyr)){install.packages("plyr"); library(plyr)}
if(!require(dplyr)){install.packages("dplyr"); library(dplyr)}
if(!require(tidyr)){install.packages("tidyr"); library(tidyr)}
if(!require(readr)){install.packages("readr"); library(readr)}
if(!require(lme4)){install.packages("lme4"); library(lme4)}
if(!require(ggplot2)){devtools::install_github("tidyverse/ggplot2"); library(ggplot2)}
if(!require(cowplot)){devtools::install_github("wilkelab/cowplot"); library(cowplot); theme_set(theme_cowplot())
  } else theme_set(theme_cowplot())
if(!require(corrplot)){install.packages("corrplot"); library(corrplot)}
if(!require(purrr)){install.packages("purrr"); library(purrr)}
if(!require(evolqg)){install.packages("evolqg"); library(evolqg)}

## Meta data

raw.mouse_meta = read_csv("./data/growth traits/F3Phenotypes_further corrected family data_corrected litter sizes.csv")
names(raw.mouse_meta) = gsub('SexAN', 'SEX', names(raw.mouse_meta))
mouse_meta = select(raw.mouse_meta, ID:COHORT)

## Markers

markers_list = llply(paste0("./data/markers/chrom", 1:19, ".csv"), read_csv)
names(markers_list) = paste0("chrom", 1:19)
for(chrom in names(markers_list)){
  names(markers_list[[chrom]])[-1] = (paste(chrom, names(markers_list[[chrom]])[-1], sep = "_"))
}

markers = Reduce(inner_join, markers_list)

chroms = seq(markers_list)
loci_per_chrom = laply(chroms, function(chrom) ncol(markers_list[[chrom]][-1])/3)

### Simulted Markers
#
#simulated_markers = llply(paste0("./data/Simulated\ independent\ genotypes/data", 1:20, ".csv"), read_csv)
#names(simulated_markers) = paste0("chrom", 1:20)
#for(chrom in names(simulated_markers)){
#  names(simulated_markers[[chrom]])[-1] = (paste(chrom, names(simulated_markers[[chrom]])[-1], sep = "_"))
#}
#
#simulated_sets = seq(simulated_markers)
#simulated_loci_per_set = laply(simulated_sets, function(chrom) ncol(simulated_markers[[chrom]][-1])/3)
#
### Simulated Chromossomes
#
#simulated_chroms_Ad = llply(paste0("./data/Simulated\ genomes/", 1:20, "_Ad.txt"), read_csv, col_names = F)
#simulated_chroms_Dm = llply(paste0("./data/Simulated\ genomes/", 1:20, "_Dom.txt"), read_csv, col_names = F)
#ID = read_csv("./data/Simulated genomes/Breed.txt", col_names = F)$X1
#for(genome in 1:20){
#  names(simulated_chroms_Ad[[genome]]) = unlist(Map(paste0, paste0(paste0("chrom", chroms), "_A"), lapply(loci_per_chrom, seq)))
#  simulated_chroms_Ad[[genome]]$ID = ID
#  names(simulated_chroms_Dm[[genome]]) = unlist(Map(paste0, paste0(paste0("chrom", chroms), "_D"), lapply(loci_per_chrom, seq)))
#  simulated_chroms_Dm[[genome]]$ID = ID
#}
#
#simulated_genomes = Map(function(x, y) inner_join(x, y, by = "ID"), 
#                        simulated_chroms_Ad, simulated_chroms_Dm)
#names(simulated_genomes[[19]])

##growth traits

raw.growth_phen = read_csv("data/growth traits/F3Phenotypes_further corrected family data_corrected litter sizes.csv")
raw.growth_phen = tbl_df(select(raw.growth_phen, c(ID, WEEK1:WEEK10)))

growth_phen = inner_join(mouse_meta, raw.growth_phen, by = "ID") %>% 
  semi_join(markers, by = "ID") %>% 
  mutate(growth12 = WEEK2 - WEEK1,
         growth23 = WEEK3 - WEEK2,
         growth34 = WEEK4 - WEEK3,
         growth45 = WEEK5 - WEEK4,
         growth56 = WEEK6 - WEEK5,
         growth67 = WEEK7 - WEEK6,
         growth78 = WEEK8 - WEEK7,
         growth89 = WEEK9 - WEEK8,
         growth910= WEEK10 - WEEK9) %>%
  select(ID, FAMILY, Dam, NURSE, SEX, LSB, LSW, COHORT, xfostpair, growth12:growth78) %>%
  na.omit %>%
  arrange(ID)

growth_markers = semi_join(markers, growth_phen, by = "ID")

growth_traits = c("growth12", "growth23", "growth34", "growth45", "growth56", "growth67", "growth78")
num_growth_traits = length(growth_traits)

m_growth_phen = gather(growth_phen, variable, value, growth12:growth78)

null.formula = "value ~ variable + variable * SEX + variable * LSB + variable * LSW + variable * COHORT"
mouse_no_fixed = lm(as.formula(null.formula), data = m_growth_phen)
m_growth_phen_std = m_growth_phen
m_growth_phen_std$value = residuals(mouse_no_fixed)

growth_phen_std = spread(m_growth_phen_std, variable, value)
growth_phen_sd = sapply(growth_phen_std[growth_traits], sd)
growth_phen_std[growth_traits] = scale(growth_phen_std[growth_traits])

# weight traits

raw.weight_phen = read_csv("data/growth traits/F3Phenotypes_further corrected family data_corrected litter sizes.csv")
raw.weight_phen = tbl_df(select(raw.weight_phen, c(ID, WEEK1:WEEK8)))

weight_phen = inner_join(mouse_meta, raw.weight_phen, by = "ID") %>% 
  semi_join(markers, by = "ID") %>% 
  select(ID, FAMILY, SEX, LSB, LSW, COHORT, WEEK1:WEEK8) %>%
  na.omit %>%
  arrange(ID)

weight_markers = semi_join(markers, weight_phen, by = "ID")

weight_traits = paste0("WEEK", 1:8)
num_weight_traits = length(weight_traits)

m_weight_phen = gather(weight_phen, variable, value, WEEK1:WEEK8)

null.formula = "value ~ variable + variable * SEX + variable * LSB + variable * LSW + variable * COHORT"
mouse_no_fixed = lm(as.formula(null.formula), data = m_weight_phen)
m_weight_phen_std = m_weight_phen
m_weight_phen_std$value = residuals(mouse_no_fixed)

weight_phen_std = spread(m_weight_phen_std, variable, value)
weight_phen_sd = sapply(weight_phen_std[weight_traits], sd)
weight_phen_std[weight_traits] = scale(weight_phen_std[weight_traits])

#Area traits

raw.area_phen = read_csv("data/area traits/areasF2F3.csv")
area_phen = inner_join(mouse_meta, raw.area_phen, by = "ID") %>% 
  semi_join(markers, by = "ID") %>% 
  select(ID, FAMILY, SEX, LSB, LSW, COHORT, area1:area7) %>%
  na.omit %>%
  arrange(ID)

area_markers = semi_join(markers, area_phen, by = "ID")

area_traits = paste0("area", 1:7)
num_area_traits = length(area_traits)

m_area_phen = gather(area_phen, variable, value, area1:area7)

null.formula = "value ~ variable * SEX + variable * LSB + variable * LSW + variable * COHORT"
mouse_no_fixed = lm(as.formula(null.formula), data = m_area_phen)
m_area_phen_std = m_area_phen
m_area_phen_std$value = residuals(mouse_no_fixed)

area_phen_std = spread(m_area_phen_std, variable, value)
area_phen_sd = sapply(area_phen_std[area_traits], sd)
area_phen_std[area_traits] = scale(area_phen_std[area_traits])

rm(list = ls(pattern='raw'))
rm(list = c("null.formula", 'mouse_no_fixed'))

##necropsy traits

raw.necropsy_phen = read_csv("data/growth traits/F3Phenotypes_further corrected family data_corrected litter sizes.csv")
raw.necropsy_phen = tbl_df(select(raw.necropsy_phen, c(ID, FATPAD, HEART:LIVER)))

necropsy_phen = inner_join(mouse_meta, raw.necropsy_phen, by = "ID") %>% 
  semi_join(markers, by = "ID") %>%
  select(ID, FAMILY, SEX, LSB, LSW, COHORT, FATPAD:LIVER) %>%
  na.omit %>%
  arrange(ID)

necropsy_markers = semi_join(markers, necropsy_phen, by = "ID")

necropsy_traits = c("FATPAD", "HEART", "KIDNEY", "LIVER", "SPLEEN")
num_necropsy_traits = length(necropsy_traits)

m_necropsy_phen = gather(necropsy_phen, variable, value, FATPAD:LIVER)

null.formula = "value ~ variable + variable * SEX + variable * LSB + variable * LSW + variable * COHORT"
mouse_no_fixed = lm(as.formula(null.formula), data = m_necropsy_phen)
m_necropsy_phen_std = m_necropsy_phen
m_necropsy_phen_std$value = residuals(mouse_no_fixed)

necropsy_phen_std = spread(m_necropsy_phen_std, variable, value)
necropsy_phen_sd = sapply(necropsy_phen_std[necropsy_traits], sd)
necropsy_phen_std[necropsy_traits] = scale(necropsy_phen_std[necropsy_traits])
