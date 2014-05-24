source('./runStanMultivariateMouseModel.R')

ggplot(m.data, aes(variable, value, color = SEX)) + geom_boxplot() + facet_wrap(~SEX)

#######################
# Classic mixed model
#######################

null.formula = "value ~ 1 + variable * SEX + variable * LSB + variable * LSW + variable * COHORT + (0 + variable|FAMILY)"
mouse.model.no.gen = lmer(as.formula(null.formula),
                          data = m.data,
                          REML = FALSE)
G_lme4 = VarCorr(mouse.model.no.gen)[['FAMILY']]
attr(G_lme4,"correlation") = NULL
dimnames(G_lme4) = list(traits, traits)

######################
# Multivariate model
######################

containsZero = function(x) ifelse(0 > x[1] & 0 < x[2], TRUE, FALSE)
find_CI = function(x, prob = 0.95){
    n = length(x)
    xs = sort(x)
    nint = floor(prob*n)
    lowest_int = abs(xs[n] - xs[1])
    for(i in 1:(n-nint)){
        current_int = abs(xs[i] - xs[i+nint])
        if(current_int <= lowest_int){
            lowest_int = current_int
            pos = i
        }
    }
    return(c(xs[pos], xs[pos+nint]))
}

mmv = stan_samples[[3]]

Gs_stan = mmv$Sigma_FAMILY
G_stan = aaply(Gs_stan, 2:3, mean)
dimnames(G_stan) = list(traits, traits)

additiveVariance = function(current_rep){
    additive_list = vector("list", num_traits)
    for(mmv in stan_samples){
        n_loci = dim(mmv$beta_addi)[2]
        addi_CI = apply(mmv$beta_addi, 2:3, find_CI)
        significant = !aaply(addi_CI, 2:3, containsZero)
        additive_matrix = mmv$beta_addi[current_rep,,]
        for(locus in 1:n_loci){
            for(trait in 1:num_traits){
                if(significant[locus, trait])
                    additive_list[[trait]] = c(additive_list[[trait]], additive_matrix[locus, trait])
            }
        }
    }
    laply(additive_list, function(x) sum(x*x)/2)
}
library(doMC)
registerDoMC(4)
aaply(1:1000, 1, additiveVariance, .parallel = TRUE)

additive_effects = adply(2:dim(mmv$beta_addi)[2], 1, function(x) data.frame(1:1000, mmv$beta_addi[,x,], locus = x))
colnames(additive_effects) = c("iterations", traits, "locus")
additive_effects = melt(additive_effects, id.vars = c('iterations','locus'))
names(additive_effects) = c('iterations', 'locus', 'trait', 'value')
ggplot(additive_effects, aes(locus, value, group = locus)) + geom_boxplot() + facet_wrap(~trait)

Gs_add = llply(stan_samples, function(mmv) aaply(mmv$beta_addi, 1, function(x) Reduce("+", alply(x, 1, function(rows) outer(rows, rows)))))
Gs = Reduce("+", Gs_add)
G_add = apply(Gs, 2:3, mean)
MatrixCompare(G_stan, G_add)
