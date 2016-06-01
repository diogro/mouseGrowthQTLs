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

raw_files = dir("data/Skull raw landmarks/")
raw_files = paste0("./data/Skull raw landmarks/", raw_files)
read_skulls = function(file){
x = scan(file, character())
char_lands = matrix(x[-c(1:3)], as.numeric(x[3]), 4, byrow = T)
lands = matrix(as.numeric(char_lands[,-1]), as.numeric(x[3]), 3)
colnames(lands) = c("x", "y", "z")
rownames(lands) = char_lands[,1]
lands}
landmarks = aperm(laply(raw_files, read_skulls), c(2, 3, 1))

library(geomorph)
raw_files[191]
skull_gpa = gpagen(landmarks)
mean_spec = findMeanSpec(skull_gpa$coords)
read.ply("./data/Multivariate/Data from Klingenberg/Ply files/1381.PLY")
