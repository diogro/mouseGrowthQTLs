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
library(ggrepel)
library(doMC)
registerDoMC(4)

source("./localShapeScripts/ggshape.R")
source("./localShapeScripts/splines.R")

raw_landmarks = read_csv("./data/Mandibles/F2F3 raw landmark positions.csv")
n_ind = dim(raw_landmarks)[1]
n_land = 15
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
for(ind in 1:length(names(n_reps))){
  landmarks_noreps[,,ind] = gpagen(landmarks[,,which(names(n_reps)[ind] == dimnames(landmarks)[[3]])])$consensus
}
dimnames(landmarks_noreps) <- list(1:15, c("x", "y"), names(n_reps))

table(raw_landmarks$ID)

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


shape = mshape
colors = rnorm(15)
ggshape_2d(mshape, wireframe, colors)
wrapMarquez(landmarks_noreps, tesselation)
