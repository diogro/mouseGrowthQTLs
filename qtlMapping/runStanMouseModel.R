library(ggplot2)
library(reshape2)
library(plyr)
library(dplyr)
library(lme4)
library(rstan)
library(gridExtra)
library(gtools)
library(glmer2stan)
library(MCMCglmm)
library(Morphometrics)

source('read.mouse.data.R')

pedigree = select(mouse.data, ID, Dam, Sire)

mouse.data = select(mouse.data, ID, FAMILY, SEX, LSB, LSW, COHORT, grow12:grow78, A1:I31)
mouse.data = mouse.data[complete.cases(mouse.data),]

num_traits = 7
traits = c( "grow12", "grow23", "grow34", "grow45", "grow56", "grow67", "grow78")

m.data = melt(mouse.data, id.vars = names(mouse.data)[c(1:6, 14:106)])
loci = as.list(m.data[grep('[AD]\\d', names(m.data))])

ggplot(m.data, aes(variable, value, color = SEX)) + geom_boxplot() + facet_wrap(~SEX)

null.formula = "value ~ 1 + variable * SEX + variable * LSB + variable * LSW + variable * COHORT + (0 + variable|FAMILY)"
mouse.model.no.gen = lmer(as.formula(null.formula),
                          data = m.data,
                          REML = FALSE)
G_lme4 = VarCorr(mouse.model.no.gen)[['FAMILY']]
attr(G_lme4,"correlation") = NULL
dimnames(G_lme4) = list(traits, traits)

null.formula = "value ~ 1 + variable * SEX + variable * LSB + variable * LSW + variable * COHORT"
mouse_no_fixed = lm(as.formula(null.formula), data = m.data)
m.data.std = m.data
m.data.std$value = residuals(mouse_no_fixed)

exclude = c(dim(m.data.std)[2]-1, dim(m.data.std)[2])
cast_formula = paste(paste(names(m.data.std[,-exclude]), collapse = " + "),
                     'variable',
                     sep = " ~ ")
mouse_data_std = dcast(m.data.std, as.formula(cast_formula))

n_locus = 1
name_locusA = paste0('A', n_locus)
name_locusD = paste0('D', n_locus)
null.formula = paste0("value ~ 1 + variable * ", name_locusA,
                      " + variable * ", name_locusD,
                      " + (variable|FAMILY)")
mouse_model = glmer2stan(as.formula(null.formula), data=m.data.std,
                         family="gaussian",
                         sample = FALSE, calcDIC = FALSE)
#write(mouse_model$model, file = "mouse_growth.stan")


N              = dim(m.data)[1]
value          = m.data.std$value
variablegrow23 = as.numeric(m.data$variable == "grow23")
variablegrow34 = as.numeric(m.data$variable == "grow34")
variablegrow45 = as.numeric(m.data$variable == "grow45")
variablegrow56 = as.numeric(m.data$variable == "grow56")
variablegrow67 = as.numeric(m.data$variable == "grow67")
variablegrow78 = as.numeric(m.data$variable == "grow78")
FAMILY         = as.integer(as.factor(m.data$FAMILY))
N_FAMILY       = length(unique(FAMILY))
mouse_list <- list(N                     = N,
                   value                 = value,
                   variablegrow23        = variablegrow23,
                   variablegrow34        = variablegrow34,
                   variablegrow45        = variablegrow45,
                   variablegrow56        = variablegrow56,
                   variablegrow67        = variablegrow67,
                   variablegrow78        = variablegrow78,
                   FAMILY                = FAMILY,
                   N_FAMILY              = N_FAMILY)
locusA = loci[[name_locusA]]
locusD = loci[[name_locusD]]
mouse_list[[name_locusA]] = locusA
mouse_list[[name_locusD]] = locusD
for(trait in grep('variablegrow', names(mouse_list), value = TRUE)){
    data_nameA = paste(trait, name_locusA, sep='_X_')
    data_nameD = paste(trait, name_locusD, sep='_X_')
    mouse_list[[ data_nameA ]] = unlist(Map(function(x, y) ifelse(x, y, 0), mouse_list[[trait]], locusA))
    mouse_list[[ data_nameD ]] = unlist(Map(function(x, y) ifelse(x, y, 0), mouse_list[[trait]], locusD))
}

mouse_stan_model = stan(file = './mouse_growth.stan', data = mouse_list, chain=1)

mm = extract(mouse_stan_model, permuted = TRUE)
extractGmat = function(x) cov(cbind(x[,1], x[,1]+x[,2], x[,1]+x[,3], x[,1]+x[,4], x[,1]+x[,5], x[,1]+x[,6], x[,1]+x[,7]))
Gs_stan = aaply(mm$vary_FAMILY, 1, extractGmat)
G_stan = aaply(Gs_stan, 2:3, mean)
dimnames(G.mcmc) = dimnames(G_stan) = list(traits, traits)
MatrixCompare(G_stan, G_lme4)
MatrixCompare(G.mcmc, G_lme4)
MatrixCompare(G.mcmc, G_stan)
cov2cor(G_stan)
cov2cor(G.mcmc)

#####################
# Multivariate model
#####################

locus_id = 1:31
name_locusA = paste0('A', locus_id)
name_locusD = paste0('D', locus_id)
N        = dim(mouse_data_std)[1]
N_locus  = length(locus_id)
K        = num_traits
mu       = numeric(K)
Sigma    = 100*diag(K)
growth   = as.matrix(mouse_data_std[grep('grow', names(mouse_data_std))])
FAMILY   = as.integer(as.factor(mouse_data_std$FAMILY))
N_FAMILY = length(unique(FAMILY))
add      = mouse_data_std[,name_locusA]
dom      = mouse_data_std[,name_locusD]
mouse_list <- list(N        = N,
                   N_locus  = N_locus,
                   K        = K,
                   mu       = mu,
                   Sigma    = Sigma,
                   growth   = growth,
                   FAMILY   = FAMILY,
                   N_FAMILY = N_FAMILY,
                   addi     = add,
                   domi     = dom)

mouse_stan_model = stan(file = './mouse_growth_multivariate_multi_locus.stan', data = mouse_list, chain=1)

mmv = extract(mouse_stan_model, permuted = TRUE)

Gs_stan = mmv$Sigma_FAMILY
G_stan = aaply(Gs_stan, 2:3, mean)
dimnames(G.mcmc) = dimnames(G_stan) = list(traits, traits)
MatrixCompare(G_stan, G_lme4)
MatrixCompare(G.mcmc, G_lme4)
MatrixCompare(G.mcmc, G_stan)
cov2cor(G_stan)
cov2cor(G.mcmc)
apply(mmv$beta_addi, 2:3, quantile, c(0.025, 0.075))
additive_effects = adply(mmv$beta_addi, 1, function(x) data.frame(x, locus = 1:31) )
additive_effects = melt(additive_effects, id.vars = c('iterations','locus'))
names(additive_effects) = c('iterations', 'locus', 'trait', 'value')
corrmat = apply(aaply(mmv$beta_addi, 1, function(x) cor(t(x))), 2:3, mean)
color2D.matplot(abs(corrmat))

ggplot(additive_effects, aes(locus, value, group = locus)) + geom_boxplot() + facet_wrap(~trait)

add_2 = (aaply(mmv$beta_addi, 1, function(x) t(x) %*% x))/2
add_variance = aaply(add_2, 2:3, mean)
