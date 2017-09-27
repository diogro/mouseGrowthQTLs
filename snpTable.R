library(dplyr)
library(readr)

snpTable = read_csv("./data/snp_table_shrinkage.csv")
snpTable %>% print(n = 219)
