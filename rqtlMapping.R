source("read_mouse_data.R")

markerPositions = as.numeric(read.csv("./data/markers/marker_positions.csv")[,3])

output = rbind(c(rep("", num_growth_traits), "", unlist(sapply(1:19, function(x) rep(x, loci_per_chrom[[x]])))),
      c(rep("", num_growth_traits), "", markerPositions),
      select(inner_join(select(growth_phen, ID, grow12:grow78, SEX), select(growth_markers, ID, matches("_A"))), -ID))
colnames(output) = c(growth_traits, "sex", names(select(growth_markers, matches("_A"))))
write.csv(output, "rqtl_input.csv", row.names = FALSE)

library(qtl)

grow = read.cross("csv", "", "rqtl_input.csv", F.gen=2, genotypes = c("-1", "0", "1"))
sex <- as.numeric(pull.pheno(grow, "sex") == "M")

out.mr <- scanone(grow, chr = 6, pheno.col = 3, method="mr")
plot(out.mr)

grow <- calc.genoprob(grow, step=1, error.prob=0.001)
out.em = scanone(grow, pheno.col = 1:7, addcovar = sex)
plot(out.em)
lodint(out.em, 6, 1.5)
