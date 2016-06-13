makeSimData = function(data, marker, percent){
  a = findA(data, marker, percent)
  marker * a
}
markerVar = function(a, data, marker){
  data = mean(data) + residuals(lm(data ~ marker))
  za = data + marker * a
  n = length(data)
  V_pa = sum((za - mean(za))^2)/(n - 1)
  gen_means = tapply(za, marker, mean)
  gen_n = table(marker) / (n - 3)
  var_marker = gen_n %*% (gen_means - mean(za))^2
  (var_marker / V_pa)[1]
}
findA = function(data, marker, percent){
  V_p = var(data)
  a = sqrt((2*percent*V_p)/(1 - percent))
  return(a)
}
makeMarker =  function(locus) data.frame(select(area_data, matches(paste0("_A", locus, "$"))))[,1]

area_data = inner_join(inner_join(area_phen, area_phen_std, 
                                  by = c("ID", "FAMILY", "SEX", "LSB", "LSW", "COHORT")) , 
                       simulated_markers[[1]][1:16], by = "ID")
sim_data = area_data %>% 
  mutate(area1 = area1.x, area6 = area6.x) %>%
  mutate(area2 = area2.x + makeSimData(marker = makeMarker(2), percent = 0.05, data = area2.y),
         area3 = area3.x + makeSimData(marker = makeMarker(2), percent = 0.05, data = area3.y),
         area4 = area4.x + makeSimData(marker = makeMarker(2), percent = 0.05, data = area4.y)) %>%
  mutate(area5 = area5.x + makeSimData(marker = makeMarker(3), percent = 0.05, data = area5.y),
         area7 = area7.x + makeSimData(marker = makeMarker(3), percent = 0.05, data = area7.y)) %>%
  mutate(area2 = area2   + makeSimData(marker = makeMarker(4), percent = 0.05, data = area2),
         area3 = area3   + makeSimData(marker = makeMarker(4), percent = 0.05, data = area3),
         area4 = area4   + makeSimData(marker = makeMarker(4), percent = 0.05, data = area4),
         area5 = area5   + makeSimData(marker = makeMarker(4), percent = 0.05, data = area5),
         area7 = area7   + makeSimData(marker = makeMarker(4), percent = 0.05, data = area7)) %>%
  select(ID:COHORT, area1:area7,
         matches(paste0("_A", 1, "$")),
         matches(paste0("_A", 2, "$")),
         matches(paste0("_A", 3, "$")),
         matches(paste0("_A", 4, "$")),
         matches(paste0("_A", 5, "$"))) %>%
  rename(A1 = sim_chrom1_A1,
         A2 = sim_chrom1_A2,
         A3 = sim_chrom1_A3,
         A4 = sim_chrom1_A4,
         A5 = sim_chrom1_A5)

only_phen = sim_data %>% select(ID:area7)
write_csv(only_phen, "~/projects/mappingTutorial/data/raw/phenotype_data.csv")

only_markers = area_data %>% 
  rename_(.dots=setNames(names(.), gsub("sim_chrom1_", "", names(.)))) %>%
  select(-(FAMILY:area7.y))
write_csv(only_markers, "~/projects/mappingTutorial/data/raw/marker_data.csv")
