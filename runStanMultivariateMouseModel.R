library(ggplot2)
library(reshape2)
library(plyr)
library(dplyr)
library(lme4)
library(rstan)
library(gridExtra)
library(gtools)
library(MCMCglmm)
library(Morphometrics)

source('read.mouse.data.R')

pedigree = select(mouse_phen, ID, Dam, Sire)

mouse_phen = select(mouse_phen, ID, FAMILY, SEX, LSB, LSW, COHORT, grow12:grow78)
complete_rows = complete.cases(mouse_phen)
mouse_phen = mouse_phen[complete_rows,]
mouse_gen = llply(mouse_gen, function(x) x[complete_rows,])

num_traits = 7
traits = c( "grow12", "grow23", "grow34", "grow45", "grow56", "grow67", "grow78")

m.data = melt(mouse_phen, id.vars = names(mouse_phen)[1:6])

null.formula = "value ~ 1 + variable * SEX + variable * LSB + variable * LSW + variable * COHORT"
mouse_no_fixed = lm(as.formula(null.formula), data = m.data)
m.data.std = m.data
m.data.std$value = residuals(mouse_no_fixed)

exclude = c(dim(m.data.std)[2]-1, dim(m.data.std)[2])
cast_formula = paste(paste(names(m.data.std[,-exclude]), collapse = " + "),
                     'variable',
                     sep = " ~ ")
mouse_phen_std = dcast(m.data.std, as.formula(cast_formula))

#####################
# Multivariate model
#####################

mouseStanInput = function(chromossome, locus_id = NULL){
    loci = grep('[AD]\\d', names(mouse_gen[[chromossome]]))
    if(is.null(locus_id))
        locus_id = 1:(length(loci)/2)
    name_locusA = paste0('A', locus_id)
    name_locusD = paste0('D', locus_id)
    N           = dim(mouse_phen_std)[1]
    N_locus     = length(locus_id)
    K           = num_traits
    mu          = numeric(K)
    Sigma       = 100*diag(K)
    growth      = as.matrix(mouse_phen_std[grep('grow', names(mouse_phen_std))])
    FAMILY      = as.integer(as.factor(mouse_phen_std$FAMILY))
    N_FAMILY    = length(unique(FAMILY))
    add         = as.matrix(mouse_gen[[chromossome]])[,name_locusA, drop = FALSE]
    dom         = as.matrix(mouse_gen[[chromossome]])[,name_locusD, drop = FALSE]
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
    return(mouse_list)
}

runStanModels = function(chromossome){
    mouse_list = mouseStanInput(chromossome)
    mouse_stan_model = stan(file = './mouse_growth_multivariate_multi_locus.stan',
                            data = mouse_list, chain=1, iter = 2000)
    return(extract(mouse_stan_model, permuted = TRUE))
}

#library(doMC)
#registerDoMC(length(mouse_gen))
#stan_samples = alply(1:length(mouse_gen), 1, runStanModels, .parallel = TRUE)
#names(stan_samples) = names(mouse_gen)
#save(stan_samples, file = "./Rdatas/stan_full_chromossome.Rdata")
#load('./Rdatas/stan_full_chromossome.Rdata')


#runStanModelsSingleLocus = function(chromossome){
    #locus_model = list()
    #n_locus = (length(mouse_gen[[chromossome]])-1)/3
    #for(i in 1:n_locus){
        #mouse_list = mouseStanInput(chromossome, locus_id = i)
        #mouse_stan_model = stan(file = './mouse_growth_multivariate_multi_locus.stan',
                                #data = mouse_list, chain=1, iter = 2000)
        #locus_model[[i]] = extract(mouse_stan_model, permuted = TRUE)
    #}
    #names(locus_model) = paste0("chrom-", chromossome, "_locus-", 1:n_locus)
    #return(locus_model)
#}
#library(doMC)
#registerDoMC(length(mouse_gen))
#stan_samples_single_locus = alply(1:length(mouse_gen), 1, runStanModelsSingleLocus, .parallel = TRUE)
#names(stan_samples_single_locus) = names(mouse_gen)
#save(stan_samples_single_locus, file = "./Rdatas/stan_single_locus_chromossome.Rdata")
load('./Rdatas/stan_single_locus_chromossome.Rdata')
