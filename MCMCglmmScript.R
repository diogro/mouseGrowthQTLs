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
                           paste(1:num.traits, 2:(num.traits+1), sep = ''), 
                           sep = ''), collapse = ', '),
               ")", sep = '')
fixed.effects = "trait:SEX - 1"

null.formula = paste(values, fixed.effects, sep = ' ~ ')

prior = list(R = list(V = diag(num.traits), n = 0.002),
             G = list(G1 = list(V = diag(num.traits) * 0.02, n = num.traits+1)))
mcmc.mouse.model = MCMCglmm( as.formula(null.formula),
                             random = ~us(trait):FAMILY,
                             data = mouse.data,
                             rcov = ~us(trait):units,
                             family = rep("gaussian", num.traits),
                             prior = prior,
                             verbose = TRUE)

G.mcmc = apply(array(mcmc.mouse.model$VCV[,1:(num.traits*num.traits)], dim = c(1000, num.traits, num.traits)), 2:3, mean)
R.mcmc = apply(array(mcmc.mouse.model$VCV[,-c(1:(num.traits*num.traits))], dim = c(1000, num.traits, num.traits)), 2:3, mean)

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

singleLocusRandomModel = function(locus, type = D, null.formula){
  mouse.data[[paste0(type, locus)]] = as.factor(mouse.data[[paste0(type, locus)]])
  
  prior = list(R = list(V = diag(num.traits), n = 0.002),
               G = list(G1 = list(V = diag(num.traits) * 0.02, n = num.traits+1),
                        G2 = list(V = diag(num.traits) * 0.02, n = num.traits+1)))
  mcmc.mouse.model = MCMCglmm( as.formula(null.formula),
                               random = ~us(at.level(D1,"1"):trait):FAMILY+
                                 us(at.level(D1,"0"):trait):FAMILY,
                               data = mouse.data,
                               rcov = ~us(trait):units,
                               family = rep("gaussian", num.traits),
                               prior = prior,
                               verbose = TRUE)
}

# all.loci.MCMC = alply(1:31, 1, runSingleLocusMCMCModel, null.formula, .progress='text')
# save(all.loci.MCMC, file= 'mouse.cromossome1.MCMC.Rdata')
load("./mouse.cromossome1.MCMC.Rdata")
