source("./read.mouse.data.R")

markers = read.csv ("./data/genotypes/marker_positions.csv", as.is = TRUE)

names_phen_markers = c(names(mouse_phen_std), markers$Marker)

N_loci = (laply(mouse_gen, length)-1)/3
chromossome_line  = c(rep("", length(names(mouse_phen_std))), unlist(Map(rep, names(mouse_gen), N_loci)))
position_line = c(rep("", length(names(mouse_phen_std))), markers[,3])

for(chrom in mouse_gen)
    for(locus in grep("A\\d", names(chrom)))
        mouse_phen_std = cbind(mouse_phen_std, chrom[[locus]])
names(mouse_phen_std) <- names_phen_markers

insertRow <- function(existingDF, newrow, r) {
  existingDF[seq(r+1,nrow(existingDF)+1),] <- existingDF[seq(r,nrow(existingDF)),]
  existingDF[r,] <- newrow
  existingDF
}

mouse_phen_gen = insertRow(mouse_phen_std, chromossome_line, 1)
mouse_phen_gen = insertRow(mouse_phen_gen, position_line, 2)
write.table(mouse_phen_gen, "./data/mouse_cross.csv", sep = ',', row.names = FALSE)

library(qtl)
library(qtlbim)
cross = read.cross("csv", dir = "data", file = "mouse_cross.csv", genotype = c(-1, 0, 1))

traitMCMC = function(trait_col){
    data_cross = qb.data(cross, pheno.col = trait_col, trait = "normal", rancov = 2)
    model_cross = qb.model(cross, epistasis = FALSE, main.nqtl = 20)
    qb.mcmc(cross, data_cross, model_cross, mydir = "mcmcqtl")
}
trait_cols = grep("grow", names(mouse_phen_gen))
#library(doMC)
#registerDoMC(4)
#trait_mcmc = llply(trait_cols, traitMCMC, .parallel = TRUE)
#save(trait_mcmc, file = "./Rdatas/qtlbim.bytrait.Rdata")
load("./Rdatas/qtlbim.bytrait.Rdata")
par(mfrow= c(2, 4))
plot(qb.scanone(trait_mcmc[[1]]))
plot(qb.scanone(trait_mcmc[[2]]))
plot(qb.scanone(trait_mcmc[[3]]))
plot(qb.scanone(trait_mcmc[[4]]))
plot(qb.scanone(trait_mcmc[[5]]))
plot(qb.scanone(trait_mcmc[[6]]))
plot(qb.scanone(trait_mcmc[[7]]))
summary(qb.scanone(trait_mcmc[[1]],type="heritability"))
