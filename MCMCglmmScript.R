library(MCMCglmm)
library(plyr)
library(reshape2)
library(lme4)
library(ggplot2)

source('read.mouse.data.R')

num.traits = 7

values = paste("cbind(", 
               paste(paste("grow", 
                           c(12, 23, 34, 45, 56, 67, 78), 
                           sep = ''), collapse = ', '),
               ")", sep = '')
fixed.effects = "trait:SEX - 1"

null.formula = paste(values, fixed.effects, sep = ' ~ ')

prior = list(R = list(V = 100*diag(7), n = 8),
             G = list(G1 = list(V = diag(7), n = 0.1)))
mcmc.mouse.model = MCMCglmm( as.formula(null.formula),
                             random = ~us(trait):FAMILY,
                             data = mouse.data,
                             rcov = ~us(trait):units,
                             family = rep("gaussian", 7),
                             prior = prior,
                             verbose = TRUE)

Gmatrix = apply(array(mcmc.mouse.model$VCV, dim = c(1000, num.traits, num.traits)), 2:3, mean)

pedigree = select(mouse.data, ID, Dam, Sire)