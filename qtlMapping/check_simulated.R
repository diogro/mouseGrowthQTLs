setwd("/home/diogro/projects/mouse-qtls")
source('read_mouse_data.R')

load("./data/Rdatas/simulated_mcmc_data2.Rdata")
load("./data/Rdatas/simulated_mcmc_data2_long.Rdata")
save(simulated_MCMC_long, file = "./data/Rdatas/simulatedEffects_data2_mcmc_long.Rdata")

effects = read_csv("./data/area\ traits/simulatedEffects_data2_mcmc.csv")
long_effects = read_csv("./data/area\ traits/simulatedEffects_data2_mcmc_long.csv")

load("./data/Rdatas/area_MCMC_null_model.Rdata")

DIC_diff = unlist(lapply(simulated_MCMC, function(x) x$DIC)) - area_MCMC_null_model$DIC
DIC_diff_long = unlist(lapply(simulated_MCMC_long, function(x) x$DIC)) - area_MCMC_null_model$DIC

ggplot(gather(data.frame(loci = 1:1000, long = DIC_diff_long, short = DIC_diff),
                         "variable", "value", long:short),
       aes(variable, value)) + geom_violin() + geom_text(aes(label = loci), position = position_jitter(width=1, height=0))

