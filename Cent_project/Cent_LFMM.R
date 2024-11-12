# lfmm2 analysis of Centaurea environmental correlation
# split into 3 groups: PNW, NE, and Europe
# use bioclim data for analysis

setwd("/gpfs1/cl/pbio3990/PopulationGenomics/metadata/")

library(LEA)
################################################################################
# Creation of a genotype matrix data file: "genotypes.lfmm"
# The data include 400 SNPs for 50 individuals.
data("tutorial")
# Write genotypes in the lfmm format
write.lfmm(tutorial.R, "genotypes.lfmm")
# Write genotypes in the geno format
write.geno(tutorial.R, "genotypes.geno")
# creation of an environment gradient file: gradient.env.
# The .env file contains a single ecological variable
# for each individual.
write.env(tutorial.C, "gradients.env")


###############################################################################
## lfmm2
###############################################################################
# load simulated data
data("offset_example")
# 200 diploid individuals genotyped at 510 SNP
Y <- offset_example$geno
# 4 environmental variables
X <- offset_example$env
mod.lfmm2 <- lfmm2(input = Y, env = X, K = 5)

# Simulate non-null effect sizes for 10 target loci
#individuals
n = 100
#loci
L = 1000
# Environmental variable
X = as.matrix(rnorm(n))
# effect sizes
B = rep(0, L)
target = sample(1:L, 10)

# GEA significance test
# showing the K = 2 estimated factors
plot(mod.lfmm2@U, col = "grey", pch = 19,
     xlab = "Factor 1",
     ylab = "Factor 2")

B[target] = runif(10, -10, 10)

# Create 3 hidden factors and their loadings
U = t(tcrossprod(as.matrix(c(-1,0.5,1.5)), X)) +
  matrix(rnorm(3*n), ncol = 3)
V <- matrix(rnorm(3*L), ncol = 3)

pv <- lfmm2.test(object = mod.lfmm2,
                 input = Y,
                 env = X,
                 full = TRUE)
plot(-log10(pv$pvalues), col = "grey", cex = .5, pch = 19)
abline(h = -log10(0.1/510), lty = 2, col = "orange")

# Simulate a matrix containing haploid genotypes
Y <- tcrossprod(as.matrix(X), B) +
  tcrossprod(U, V) +
  matrix(rnorm(n*L, sd = .5), nrow = n)
Y <- matrix(as.numeric(Y > 0), ncol = L)

# Fitting an LFMM with K = 3 factors
mod <- lfmm2(input = Y, env = X, K = 3)

# Computing P-values and plotting their minus log10 values
pv <- lfmm2.test(object = mod,
                 input = Y,
                 env = X,
                 linear = TRUE)
plot(-log10(pv$pvalues), col = "grey", cex = .6, pch = 19)
points(target, -log10(pv$pvalues[target]), col = "red")


