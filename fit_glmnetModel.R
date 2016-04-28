source("./read.mouse.data.R")

if(!require(glmnet)) {install.packages("glmnet"); library(glmnet)}

mouse_phen_std
