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

Gs_stan = mmv$Sigma_FAMILY
G_stan = aaply(Gs_stan, 2:3, mean)
dimnames(G_stan) = list(traits, traits)
apply(mmv$beta_addi, 2:3, quantile, c(0.025, 0.975))
#additive_effects = adply(mmv$beta_addi, 2:3, function(x) data.frame(x, locus = 10) )
additive_effects = data.frame(1:1000, mmv$beta_addi[,1,], locus = 10)
colnames(additive_effects) = c("iterations", paste0("trait", 1:7), "locus")
additive_effects = melt(additive_effects, id.vars = c('iterations','locus'))
names(additive_effects) = c('iterations', 'locus', 'trait', 'value')
corrmat = apply(aaply(mmv$beta_addi, 1, function(x) cor(t(x))), 2:3, mean)
color2D.matplot(abs(corrmat))

ggplot(additive_effects, aes(trait, value, group = trait)) + geom_boxplot()# + facet_wrap(~trait)

mean_additive = aaply(mmv$beta_addi, 2:3, mean)
add_2 = outer(mean_additive, mean_additive)/2
