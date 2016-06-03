library(evolqg)
library(plyr)
library(doMC)
registerDoMC(4)

# source("./ggshape.R")
# source("./splines.R")
#
# load("./Aux.RData")
# load("~/Dropbox/labbio/data/cov_bayes_data/Rdatas/monkeys.RData")
#
# OTU = main.data$Homo_sapiens
# tesselation = Aux$single.tessel.38
#
# wrapMarquez(OTU$ss, tesselation)

raw_files = dir("data/Skull raw landmarks/")
raw_files = paste0("./data/Skull raw landmarks/", raw_files)
read_skulls = function(file){
  x = scan(file, character())
  char_lands = matrix(x[-c(1:3)], as.numeric(x[3]), 4, byrow = T)
  lands = matrix(as.numeric(char_lands[,-1]), as.numeric(x[3]), 3)
  colnames(lands) = c("x", "y", "z")
  rownames(lands) = char_lands[,1]
  lands
}
raw_landmarks = laply(raw_files, read_skulls)
landmarks = aperm(raw_landmarks, c(2, 3, 1))

library(geomorph)

skull_gpa = gpagen(landmarks)
ref = mshape(skull_gpa$coords)
mean_spec = findMeanSpec(skull_gpa$coords)
raw_files[[mean_spec]]
raw_mesh = read.ply("./data/landmarks/1381.PLY")
mean_mesh = warpRefMesh(raw_mesh, landmarks[,,mean_spec], ref)

pc.res <- prcomp(two.d.array(skull_gpa$coords))
pcdata <- pc.res$x
rotation <- pc.res$rotation
k <- dim(skull_gpa$coords)[2] ; p <- dim(skull_gpa$coords)[1]
pcaxis1.quatiles <- quantile(pcdata[, 1], seq(0, 1, length.out = 10))
pcaxis2.quatiles <- quantile(pcdata[, 2], seq(0, 1, length.out = 10))

pc1_loading = pcaxis1.quatiles[[10]]
pc2_loading = pcaxis2.quatiles[[10]]

pc1_shape <- function(i) {
  pc.vec <- rep(0, dim(pcdata)[2])
  pc.vec[1] <- pcaxis1.quatiles[[i]]
  pc.vec[2] <- pcaxis2.quatiles[[i]]
  target <- arrayspecs(as.matrix(pc.vec %*% (t(rotation))), p,k)[,,1] + ref
  warpRefMesh(raw_mesh, landmarks[,,mean_spec], target)
}
pc1_shapes_list = llply(1:10, pc1_shape, .progress = "text")
shape = pc1_shapes_list[[1]]
plot_shape = function(shape){
  open3d()
  shade3d(rotate3d(shape, 8/7*pi, 1, 0, 0))
}

for(i in 1:10){
plot_shape(pc1_shapes_list[[i]])
rgl.postscript(paste0("./animation/pc1_animation_", i, ".eps"), fmt="eps")
snapshot3d(paste0("./animation/pc1_animation_", i, ".png"))
}

covm = cov(two.d.array(skull_gpa$coords))
MeanMatrixStatistics(covm)
