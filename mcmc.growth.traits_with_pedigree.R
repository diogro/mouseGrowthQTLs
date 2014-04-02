library(MCMCglmm)

source('~/projects/mouse-qtls/read.mouse.data.R')

ped<-data.frame(read.csv(file="~/projects/mouse-qtls/data/F3_pedigree.csv", header=TRUE))

names(data) <- gsub('SexAN', 'SEX', names(data))

n.t <- 7

prior <- list(R = list(V = diag(7), n = 0.002),
             G = list(G1 = list(V = diag(7) * 0.02, n = 8)))

mcmc.growth.traits <- MCMCglmm(cbind(grow12, grow23, grow34, grow45, grow56, grow67, grow78) ~ trait:SEX - 1,
                       random = ~us(trait):animal,
                       data = data,
                       rcov = ~us(trait):units,
                       family = rep("gaussian", n.t),
                       pedigree = ped,
                       nitt=1300000, thin = 1000, burnin = 300000,
                       prior = prior,
                       verbose = TRUE)
