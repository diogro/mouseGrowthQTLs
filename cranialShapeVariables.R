library(evolqg)
library(doMC)
registerDoMC(4)

source("./ggshape.R")
source("./splines.R")

load("./Aux.RData")
load("~/Dropbox/labbio/data/cov_bayes_data/Rdatas/monkeys.RData")

OTU = main.data$Homo_sapiens
tesselation = Aux$single.tessel.38

wrapMarquez(OTU$ss, tesselation)

