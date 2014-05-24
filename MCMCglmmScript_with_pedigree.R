library(MCMCglmm)
library(Morphometrics)

source("./MCMCglmmScript.R")
source('./read.mouse.data.R')

pedigree <- data.frame(read.csv(file="./data/F3_pedigree.csv", header=TRUE))
pedigree = pedigree[-c(1997:2000),]

mouse.data = select(mouse.data, ID, FAMILY, SEX, LSB, LSW, COHORT, grow12:grow78, A1:I31)
mouse.data = mouse.data[complete.cases(mouse.data[,c('ID', 'FAMILY', 'SEX', 'LSB', 'LSW', 'COHORT')]),]
mouse.data$animal = mouse.data$ID
mouse.data = mouse.data[mouse.data$animal %in% pedigree$ID,]
#pedigree = pedigree[pedigree$ID %in% mouse.data$animal,]
num.traits = 7

sum((na.omit(pedigree[, 2]) %in% pedigree[, 1]) == FALSE) > 0 & any(is.na(pedigree[, 2]))
sum((na.omit(pedigree[, 3]) %in% pedigree[, 1]) == FALSE) > 0 & any(is.na(pedigree[, 3]))

value = paste("cbind(",
               paste(paste("grow",
                           paste(1:num.traits, 2:(num.traits+1), sep = ''),
                           sep = ''), collapse = ', '),
               ")", sep = '')

fixed.effects = "trait:SEX + trait:LSB + trait:LSW + trait:COHORT - 1"

null.formula = paste(value, fixed.effects, sep = ' ~ ')

mcmc.mouse.model = runNullMCMCModel(null.formula, pl=TRUE)
summary(mcmc.mouse.model)

G.mcmc = apply(array(mcmc.mouse.model$VCV[,1:(num.traits*num.traits)], dim = c(1000, num.traits, num.traits)), 2:3, mean)
R.mcmc = apply(array(mcmc.mouse.model$VCV[,-c(1:(num.traits*num.traits))], dim = c(1000, num.traits, num.traits)), 2:3, mean)

prior <- list(R = list(V = diag(7), n = 0.002),
              G = list(G1 = list(V = diag(7) * 0.02, n = 8)))

start <- list(R = list(V = R.mcmc), G = list(G1 = G.mcmc))#, liab = matrix(mcmc.mouse.model$Liab[1,], ncol = num.traits))

mcmc.growth.pedigree <- MCMCglmm(as.formula(null.formula),
                                 random = ~us(trait):animal,
                                 data = mouse.data,
                                 rcov = ~us(trait):units,
                                 family = rep("gaussian", num.traits),
                                 pedigree = pedigree,
                                 start = t(start),
                                 nitt=11000, thin = 10, burnin = 1000,
                                 prior = prior,
                                 verbose = TRUE)
summary(mcmc.growth.pedigree)
G.mcmc.ped = apply(array(mcmc.growth.pedigree$VCV[,1:(num.traits*num.traits)], dim = c(1000, num.traits, num.traits)), 2:3, mean)
R.mcmc.ped = apply(array(mcmc.growth.pedigree$VCV[,-c(1:(num.traits*num.traits))], dim = c(1000, num.traits, num.traits)), 2:3, mean)
round(cov2cor(G.mcmc.ped), 3)
round(cov2cor(G.mcmc), 3)
plot(1:10)

MatrixCompare(G.mcmc.ped, G.mcmc)

ids = as.list(mouse.data$ID)
pedigree_family = list()
while(length(ids) > 0){
    mouse = ids[[1]]
    dam = filter(pedigree, ID == mouse)$Dam
    family = filter(pedigree, Dam == dam)$ID
    pedigree_family[[length(pedigree_family)+1]] = family
    ids[ids %in% family] = NULL
}
mouse_family = dlply(mouse.data, 'FAMILY', function(x) x$ID)
mouse_family = llply(mouse_family, sort)
pedigree_family = llply(pedigree_family, sort)
pedigree_family = pedigree_family[order(laply(pedigree_family, '[', 1))]
mouse_family = mouse_family[order(laply(mouse_family, '[', 1))]
all.equal(mouse_family, pedigree_family)
