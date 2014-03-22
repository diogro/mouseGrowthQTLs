library(MCMCglmm)
library(plyr)
library(reshape2)
library(lme4)
library(ggplot2)

source('read.mouse.data.R')

pedigree = select(mouse.data, ID, Dam, Sire)

num.traits = 7

values = paste("cbind(", 
               paste(paste("grow", 
                           paste(1:7, 2:8, sep = ''), 
                           sep = ''), collapse = ', '),
               ")", sep = '')
fixed.effects = "trait:SEX - 1"

null.formula = paste(values, fixed.effects, sep = ' ~ ')

prior = list(R = list(V = diag(7), n = 0.002),
             G = list(G1 = list(V = diag(7) * 0.02, n = 8)))
mcmc.mouse.model = MCMCglmm( as.formula(null.formula),
                             random = ~us(trait):FAMILY,
                             data = mouse.data,
                             rcov = ~us(trait):units,
                             family = rep("gaussian", 7),
                             prior = prior,
                             verbose = TRUE)

G.mcmc = apply(array(mcmc.mouse.model$VCV, dim = c(1000, num.traits, num.traits)), 2:3, mean)

runSingleLocusMCMCModel <- function(locus, null.formula){
  genotype.formula = paste(null.formula, 
                           paste(paste('trait:', c('A', 'D', 'I'),
                                       locus, sep = ''), collapse = ' + '),
                           sep = ' + ')
  prior = list(R = list(V = diag(7), n = 0.002),
               G = list(G1 = list(V = diag(7) * 0.02, n = 8)))
  mcmc.mouse.model = MCMCglmm( as.formula(genotype.formula),
                               random = ~us(trait):FAMILY,
                               data = mouse.data,
                               rcov = ~us(trait):units,
                               family = rep("gaussian", 7),
                               prior = prior,
                               verbose = FALSE)
  return(mcmc.mouse.model)
}

# all.loci.MCMC = alply(1:31, 1, runSingleLocusMCMCModel, null.formula, .progress='text')
# save(all.loci.MCMC, file= 'mouse.cromossome1.MCMC.Rdata')
load("./mouse.cromossome1.MCMC.Rdata")