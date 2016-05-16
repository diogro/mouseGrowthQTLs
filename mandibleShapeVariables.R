library(evolqg)
library(readr)
library(plyr)
library(dplyr)
library(reshape2)
library(tidyr)
library(geomorph)
library(expm)
library(numDeriv)
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(ggrepel)
library(doMC)
registerDoMC(4)
if(!require(wesanderson)) { install.packages("wesanderson"); library(wesanderson) }


source("./localShapeScripts/ggshape.R")
source("./localShapeScripts/splines.R")

raw_landmarks = read_csv("./data/Mandibles/F2F3 raw landmark positions.csv")
n_ind = dim(raw_landmarks)[1]
n_land = 15
n_reps = table(raw_landmarks$ID)

landmarks = array(NA, dim = c(n_land, 2, n_ind))
for(i in seq(n_ind)){
  for(j in seq(n_land))
    landmarks[j,,i] = as.numeric(raw_landmarks[i, paste0(c("X", "Y"), j)])
}
dimnames(landmarks) <- list(1:15, c("x", "y"), raw_landmarks$ID)

for(ind in names(n_reps)){
  for(i in which(ind == dimnames(landmarks)[[3]])){
    if(all(landmarks[1,1,i] > landmarks[-1,1,i])) landmarks[,1,i] = -landmarks[,1,i]
  }
}

pdf("mandible_plots.pdf", onefile = TRUE)
for(ind in names(n_reps)){
  par(mfrow = c(2, n_reps[ind]/2))
  for(i in which(ind == dimnames(landmarks)[[3]])){
    plot(landmarks[,,i], main = ind)
    lines(landmarks[,,i])
  }
}
dev.off()

landmarks_noreps = array(NA, dim = c(n_land, 2, length(names(n_reps))))
size = array(NA, dim = c(length(names(n_reps))))
for(ind in 1:length(names(n_reps))){
  gpa = gpagen(landmarks[,,which(names(n_reps)[ind] == dimnames(landmarks)[[3]])])
  landmarks_noreps[,,ind] = gpa$consensus
  size[ind] = mean(gpa$Csize)
}
dimnames(landmarks_noreps) <- list(1:15, c("x", "y"), names(n_reps))
names(size) = names(n_reps)

gpa = gpagen(landmarks_noreps)
gpa$coords
mshape = mshape(gpa$coords)
plot(gpa)

tesselation = matrix(
  c(1, 2, 15,
    2, 15, 14,
    2, 14, 13, 
    2, 3, 13,
    3, 12, 13,
    3, 4, 12,
    4, 5, 6,
    4, 6, 12,
    6, 9, 12, 
    6, 7, 9,
    7, 8, 9,
    9, 10, 11,
    9, 11, 12), ncol = 3, byrow = TRUE)

wireframe = matrix(
  c(1, 2,
    1, 15, 
    2, 15,
    2, 14,
    14, 15,
    14, 13,
    2, 3,
    3, 13,
    2, 13,
    3, 4,
    3, 12,
    4, 12,
    4, 6,
    12, 13,
    6, 12,
    5, 6,
    4, 5,
    6, 7,
    6, 9,
    7, 9,
    7, 8,
    8, 9,
    9, 10,
    10, 11,
    9, 11,
    9, 12,
    11, 12), ncol = 2, byrow = TRUE)



#shape_variables = wrapMarquez(landmarks_noreps, tesselation)
#save(shape_variables, file = "./data/Rdatas/mandibleShape.Rdata")
load("./data/Rdatas/mandibleShape.Rdata")
shape_variables$cs = size


shape_cor = cor(shape_variables$local)
eigen_1 = eigen(shape_cor)$vectors[,1]
plotMatrix(shape_cor)

colors = rnorm(15)
positions = data.frame(
  mshape[t(tesselation),], 
  id = rep(1:13, each = 3))

centroids = ddply(positions, .(id), numcolwise(mean)) 

values = data.frame(id = 1:13, eigen_1 = shape_variables$local[16,])
datapoly <- merge(values, positions, by=c("id"))

min_s = min(shape_variables$local)
max_s = max(shape_variables$local)

min_size = min(shape_variables$cs)
max_size = max(shape_variables$cs)

mypalette = colorRampPalette(c("blue", "white", "red"))
(p <- ggplot(datapoly, aes(x=x, y=y)) + 
  geom_polygon(color = "black", aes(fill=eigen_1, group=id)) + 
  geom_text(data = centroids, aes(x, y, label = id)) +
  scale_fill_gradientn('Shape\nVariables', colours = mypalette(10), limits=c(min_s, max_s)) + 
  scale_alpha_continuous(guide = FALSE) + 
  theme_shape())
save_plot("mandible_loadings.png", p, base_height = 5, base_aspect_ratio = 2)

ggshape_2d(shape_variables$reference, wireframe, colors)
