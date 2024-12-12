# Required packages
# Adjusting genotype-environment association models
library(LEA)
library(tidyverse)
library(vcfR)
library(SNPfiltR)

library(conflicted)  

library(tidyverse)
conflict_prefer("filter", "dplyr")
conflict_prefer("lag", "dplyr")

options(bitmapType = "cairo")

setwd("~/projects/eco_genomics/population_genomics/")

vcf <- read.vcfR("~/projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.vcf.gz")

meta <- read.csv("/gpfs1/cl/pbio3990/GroupProjects/Cent_climadapt/metaA.csv", row.names = "X")
dim(meta)

meta2 <- meta[meta$id %in% colnames(vcf@gt[, -1]),]
dim(meta2)

vcf <- vcf[,colnames(vcf@gt[,-1]) %in% meta$id]
missing_by_snp(vcf)
vcf <-missing_by_snp(vcf, 0.75)

#write.vcf(vcf,"~/projects/eco_genomics/Cent_project/vcf_final.filtered_setA.vcf.gz")
#system("gunzip -c ~/projects/eco_genomics/Cent_project/vcf_final.filtered_setA.vcf.gz > ~/vcf_final.filtered_setA.vcf")

#vcf2lfmm("/gpfs1/home/l/e/lericsso/vcf_final.filtered_setA.vcf","/gpfs1/home/l/e/lericsso/projects/eco_genomics/Cent_project/vcf_final.filtered_setA.lfmm")
vcflfmm <- read.lfmm("/users/l/e/lericsso/vcf_final.filtered_SRK2.lfmm")
#Unnecessary if using vcf2lfmm:
#geno <- vcf2geno(input.file = "/gpfs1/home/l/e/lericsso/vcf_final.filtered_setA.vcf", output.file = "outputs/vcf_final.filtered_setA.geno")
#geno_vcf <- read.vcfR("/users/l/e/lericsso/projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.vcf.gz")

#######
#climate data
#######
#unnecessary if using vcf2lfmm
#geno_lfmm <- geno2lfmm(input.file = "/users/l/e/lericsso/projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.geno", output.file = "/users/l/e/lericsso/projects/eco_genomics/Cent_project/geno.lfmm", force = TRUE)

# All bioclim variables
write.env(meta2[,9:27], output.file = "/users/l/e/lericsso/projects/eco_genomics/Cent_project/climate.env")
meta_env <- read.env("/users/l/e/lericsso/projects/eco_genomics/Cent_project/climate.env")

# 1. TempM
write.env(meta2[,9], output.file = "/users/l/e/lericsso/projects/eco_genomics/Cent_project/TempM.env")
TempM_env <- read.env("/users/l/e/lericsso/projects/eco_genomics/Cent_project/TempM.env")

# 2. TempR
write.env(meta2[,15], output.file = "/users/l/e/lericsso/projects/eco_genomics/Cent_project/TempR.env")
TempR_env <- read.env("/users/l/e/lericsso/projects/eco_genomics/Cent_project/TempR.env")


# 3. Prec
write.env(meta2[,20], output.file = "/users/l/e/lericsso/projects/eco_genomics/Cent_project/Prec.env")
Prec_env <- read.env("/users/l/e/lericsso/projects/eco_genomics/Cent_project/Prec.env")

# 4. MDiurnalR
write.env(meta2[,10], output.file = "/users/l/e/lericsso/projects/eco_genomics/Cent_project/MDiurnalR.env")
MDiurnalR_env <- read.env("/users/l/e/lericsso/projects/eco_genomics/Cent_project/MDiurnalR.env")

# 5. PrecDM
write.env(meta2[,22], output.file = "/users/l/e/lericsso/projects/eco_genomics/Cent_project/PrecDM.env")
PrecDM_env <- read.env("/users/l/e/lericsso/projects/eco_genomics/Cent_project/PrecDM.env")

# 6. PrecWmQ
write.env(meta2[,26], output.file = "/users/l/e/lericsso/projects/eco_genomics/Cent_project/PrecWmQ.env")
PrecWmQ_env <- read.env("/users/l/e/lericsso/projects/eco_genomics/Cent_project/PrecWmQ.env")

########
# Admixture, population differentiation tests, missing genotype imputation
########
# main options
# K = number of ancestral populations
# entropy = TRUE computes the cross-entropy criterion,
# CPU = 4 is the number of CPU used (hidden input)
project = NULL
project = snmf("/users/l/e/lericsso/vcf_final.filtered_SRK2.lfmm",
               K = 1:10,
               entropy = TRUE,
               repetitions = 10,
               project = "new")

plot(project, col = "blue", pch = 19, cex = 1.2)
best = which.min(cross.entropy(project, K = 5))
########
#color palette
########
library(RColorBrewer)
cols <- brewer.pal(6, "Set3")
cols
pal <- colorRampPalette(cols)

palette(brewer.pal(n = 6, name = "Set2"))

######
p = snmf.pvalues(project,
                 entropy = TRUE,
                 ploidy = 2,
                 K = 5)
pvalues = p$pvalues
par(mfrow = c(2,1))
hist(pvalues, col = "orange")
plot(-log10(pvalues), pch = 19, col = "blue", cex = .5)
dev.off()
#######

project.missing = snmf("/users/l/e/lericsso/vcf_final.filtered_SRK2.lfmm", K = 5,
                       entropy = TRUE, repetitions = 10,
                       project = "new")

# select the run with the lowest cross-entropy value
best = which.min(cross.entropy(project.missing, K = 5))
# Impute the missing genotypes
LEA::impute(project.missing, "/users/l/e/lericsso/vcf_final.filtered_SRK2.lfmm",
            method = 'mode', K = 5, run = best)
## Missing genotype imputation for K = 5
## Missing genotype imputation for run = 5
## Results are written in the file: /users/l/e/lericsso/projects/eco_genomics/Cent_project/geno.lfmm_imputed.lfmm
# Proportion of correct imputation results
#dat.imp = read.lfmm("/users/l/e/lericsso/projects/eco_genomics/Cent_project/geno.lfmm_imputed.lfmm")
#mean("geno.lfmm"["geno.lfmm" == 9] == dat.imp["geno.lfmm_imputed.lfmm" == 9] )

dat.imp = read.lfmm("/users/l/e/lericsso/vcf_final.filtered_SRK2.lfmm_imputed.lfmm")

######
# LFMM2 Analysis Mean Temp
######

Y <- dat.imp
dim(Y)
setAloci <- which(colMeans(Y)!=0 & colMeans(Y)!=2)
Y <- dat.imp[,setAloci]
dim(Y)
X <- meta2$TempM


mod.lfmm2 <- lfmm2(Y, X, K = 5)


# Simulate non-null effect sizes for 10 target loci
#individuals
dim(meta2)
n = 118
#loci
L = 11238

# Environmental variable
X = as.matrix(rnorm(n))
# effect sizes
B = rep(0, L)
target = sample(1:L, 10)

B[target] = runif(10, -10, 10)
# Create 3 hidden factors and their loadings
U = t(tcrossprod(as.matrix(c(-1,0.5,1.5)), X)) +
  matrix(rnorm(3*n), ncol = 3)
V <- matrix(rnorm(3*L), ncol = 3)

pv <- lfmm2.test(mod.lfmm2,
                 Y,
                 X,
                 full = TRUE)
plot(-log10(pv$pvalues), col = "grey", cex = .5, pch = 19)
abline(h = -log10(0.1/510), lty = 2, col = "orange")

# GEA significance test
# showing the K = 5 estimated factors
plot(mod.lfmm2@U, pch = 19,
     xlab = "Latent Factors",
     ylab = "Mean Temperature",
     col = as.factor(meta2$region))
legend("bottomleft", legend= levels(as.factor(meta2$region)), pch=16, col= unique(as.factor(meta2$region)))

# Simulate a matrix containing haploid genotypes
Y <- tcrossprod(as.matrix(X), B) +
  tcrossprod(U, V) +
  matrix(rnorm(n*L, sd = .5), nrow = n)
Y <- matrix(as.numeric(Y > 0), ncol = L)

# Fitting an LFMM with K = 5 factors
mod <- lfmm2(Y, TempM_env, K = 5)

# Computing P-values and plotting their minus log10 values
pv <- lfmm2.test(mod,
                 Y,
                 TempM_env,
                 linear = TRUE)
plot(-log10(pv$pvalues), col = "grey", cex = .6, pch = 19)
points(target, -log10(pv$pvalues[target]), col = "red")

######
# LFMM2 Analysis Temp Range
######

Y <- dat.imp
dim(Y)
setAloci <- which(colMeans(Y)!=0 & colMeans(Y)!=2)
Y <- dat.imp[,setAloci]
dim(Y)
X <- meta2$TempR


mod.lfmm2 <- lfmm2(Y, X, K = 5)


# Simulate non-null effect sizes for 10 target loci
#individuals
dim(meta2)
n = 118
#loci
L = 11238

# Environmental variable
X = as.matrix(rnorm(n))
# effect sizes
B = rep(0, L)
target = sample(1:L, 10)

B[target] = runif(10, -10, 10)
# Create 3 hidden factors and their loadings
U = t(tcrossprod(as.matrix(c(-1,0.5,1.5)), X)) +
  matrix(rnorm(3*n), ncol = 3)
V <- matrix(rnorm(3*L), ncol = 3)

pv <- lfmm2.test(mod.lfmm2,
                 Y,
                 X,
                 full = TRUE)
plot(-log10(pv$pvalues), col = "grey", cex = .5, pch = 19)
abline(h = -log10(0.1/510), lty = 2, col = "orange")

# GEA significance test
# showing the K = 5 estimated factors
plot(mod.lfmm2@U, pch = 19,
     xlab = "Latent Factors",
     ylab = "Temperature Range",
     col = as.factor(meta2$region))
legend("bottomleft", legend= levels(as.factor(meta2$region)), pch=16, col=unique(as.factor(meta2$region)))

# Simulate a matrix containing haploid genotypes
Y <- tcrossprod(as.matrix(X), B) +
  tcrossprod(U, V) +
  matrix(rnorm(n*L, sd = .5), nrow = n)
Y <- matrix(as.numeric(Y > 0), ncol = L)

# Fitting an LFMM with K = 5 factors
mod <- lfmm2(Y, TempR_env, K = 5)

# Computing P-values and plotting their minus log10 values
pv <- lfmm2.test(mod,
                 Y,
                 TempR_env,
                 linear = TRUE)
plot(-log10(pv$pvalues), col = "grey", cex = .6, pch = 19)
points(target, -log10(pv$pvalues[target]), col = "red")

######
# LFMM2 Analysis Prec
######

Y <- dat.imp
dim(Y)
setAloci <- which(colMeans(Y)!=0 & colMeans(Y)!=2)
Y <- dat.imp[,setAloci]
dim(Y)
X <- meta2$Prec


mod.lfmm2 <- lfmm2(Y, X, K = 5)


# Simulate non-null effect sizes for 10 target loci
#individuals
dim(meta2)
n = 118
#loci
L = 11238

# Environmental variable
X = as.matrix(rnorm(n))
# effect sizes
B = rep(0, L)
target = sample(1:L, 10)

B[target] = runif(10, -10, 10)
# Create 3 hidden factors and their loadings
U = t(tcrossprod(as.matrix(c(-1,0.5,1.5)), X)) +
  matrix(rnorm(3*n), ncol = 3)
V <- matrix(rnorm(3*L), ncol = 3)

pv <- lfmm2.test(mod.lfmm2,
                 Y,
                 X,
                 full = TRUE)
plot(-log10(pv$pvalues), col = "grey", cex = .5, pch = 19)
abline(h = -log10(0.1/510), lty = 2, col = "orange")

# GEA significance test
# showing the K = 5 estimated factors
plot(mod.lfmm2@U, pch = 19,
     xlab = "Latent Factors",
     ylab = "Precipitation",
     col = as.factor(meta2$region))
legend("bottomleft", legend= levels(as.factor(meta2$region)), pch=16, col=unique(as.factor(meta2$region)))

# Simulate a matrix containing haploid genotypes
Y <- tcrossprod(as.matrix(X), B) +
  tcrossprod(U, V) +
  matrix(rnorm(n*L, sd = .5), nrow = n)
Y <- matrix(as.numeric(Y > 0), ncol = L)

# Fitting an LFMM with K = 5 factors
mod <- lfmm2(Y, Prec_env, K = 5)

# Computing P-values and plotting their minus log10 values
pv <- lfmm2.test(mod,
                 Y,
                 Prec_env,
                 linear = TRUE)
plot(-log10(pv$pvalues), col = "grey", cex = .6, pch = 19)
points(target, -log10(pv$pvalues[target]), col = "red")

######
# LFMM2 Analysis MDiurnalR
######

Y <- dat.imp
dim(Y)
setAloci <- which(colMeans(Y)!=0 & colMeans(Y)!=2)
Y <- dat.imp[,setAloci]
dim(Y)
X <- meta2$MDiurnalR


mod.lfmm2 <- lfmm2(Y, X, K = 5)


# Simulate non-null effect sizes for 10 target loci
#individuals
dim(meta2)
n = 118
#loci
L = 11238

# Environmental variable
X = as.matrix(rnorm(n))
# effect sizes
B = rep(0, L)
target = sample(1:L, 10)

B[target] = runif(10, -10, 10)
# Create 3 hidden factors and their loadings
U = t(tcrossprod(as.matrix(c(-1,0.5,1.5)), X)) +
  matrix(rnorm(3*n), ncol = 3)
V <- matrix(rnorm(3*L), ncol = 3)

pv <- lfmm2.test(mod.lfmm2,
                 Y,
                 X,
                 full = TRUE)
plot(-log10(pv$pvalues), col = "grey", cex = .5, pch = 19)
abline(h = -log10(0.1/510), lty = 2, col = "orange")

# GEA significance test
# showing the K = 5 estimated factors
plot(mod.lfmm2@U, pch = 19,
     xlab = "Latent Factors",
     ylab = "Mean Diurnal Range",
     col = as.factor(meta2$region))
legend("bottomleft", legend= levels(as.factor(meta2$region)), pch=16, col=unique(as.factor(meta2$region)))

# Simulate a matrix containing haploid genotypes
Y <- tcrossprod(as.matrix(X), B) +
  tcrossprod(U, V) +
  matrix(rnorm(n*L, sd = .5), nrow = n)
Y <- matrix(as.numeric(Y > 0), ncol = L)

# Fitting an LFMM with K = 5 factors
mod <- lfmm2(Y, MDiurnalR_env, K = 5)

# Computing P-values and plotting their minus log10 values
pv <- lfmm2.test(mod,
                 Y,
                 MDiurnalR_env,
                 linear = TRUE)
plot(-log10(pv$pvalues), col = "grey", cex = .6, pch = 19)
points(target, -log10(pv$pvalues[target]), col = "red")

######
# LFMM2 Analysis PrecDM
######

Y <- dat.imp
dim(Y)
setAloci <- which(colMeans(Y)!=0 & colMeans(Y)!=2)
Y <- dat.imp[,setAloci]
dim(Y)
X <- meta2$PrecDM


mod.lfmm2 <- lfmm2(Y, X, K = 5)


# Simulate non-null effect sizes for 10 target loci
#individuals
dim(meta2)
n = 118
#loci
L = 11238

# Environmental variable
X = as.matrix(rnorm(n))
# effect sizes
B = rep(0, L)
target = sample(1:L, 10)

B[target] = runif(10, -10, 10)
# Create 3 hidden factors and their loadings
U = t(tcrossprod(as.matrix(c(-1,0.5,1.5)), X)) +
  matrix(rnorm(3*n), ncol = 3)
V <- matrix(rnorm(3*L), ncol = 3)

pv <- lfmm2.test(mod.lfmm2,
                 Y,
                 X,
                 full = TRUE)
plot(-log10(pv$pvalues), col = "grey", cex = .5, pch = 19)
abline(h = -log10(0.1/510), lty = 2, col = "orange")

# GEA significance test
# showing the K = 5 estimated factors
plot(mod.lfmm2@U, pch = 19,
     xlab = "Latent Factors",
     ylab = "Precipitation During Driest Month",
     col = as.factor(meta2$region))
legend("bottomleft", legend= levels(as.factor(meta2$region)), pch=16, col=unique(as.factor(meta2$region)))

# Simulate a matrix containing haploid genotypes
Y <- tcrossprod(as.matrix(X), B) +
  tcrossprod(U, V) +
  matrix(rnorm(n*L, sd = .5), nrow = n)
Y <- matrix(as.numeric(Y > 0), ncol = L)

# Fitting an LFMM with K = 5 factors
mod <- lfmm2(Y, PrecDM_env, K = 5)

# Computing P-values and plotting their minus log10 values
pv <- lfmm2.test(mod,
                 Y,
                 PrecDM_env,
                 linear = TRUE)
plot(-log10(pv$pvalues), col = "grey", cex = .6, pch = 19)
points(target, -log10(pv$pvalues[target]), col = "red")

######
# LFMM2 Analysis PrecWmQ
######

Y <- dat.imp
dim(Y)
setAloci <- which(colMeans(Y)!=0 & colMeans(Y)!=2)
Y <- dat.imp[,setAloci]
dim(Y)
X <- meta2$PrecWmQ


mod.lfmm2 <- lfmm2(Y, X, K = 5)


# Simulate non-null effect sizes for 10 target loci
#individuals
dim(meta2)
n = 118
#loci
L = 11238

# Environmental variable
X = as.matrix(rnorm(n))
# effect sizes
B = rep(0, L)
target = sample(1:L, 10)

B[target] = runif(10, -10, 10)
# Create 3 hidden factors and their loadings
U = t(tcrossprod(as.matrix(c(-1,0.5,1.5)), X)) +
  matrix(rnorm(3*n), ncol = 3)
V <- matrix(rnorm(3*L), ncol = 3)

pv <- lfmm2.test(mod.lfmm2,
                 Y,
                 X,
                 full = TRUE)
plot(-log10(pv$pvalues), col = "grey", cex = .5, pch = 19)
abline(h = -log10(0.1/510), lty = 2, col = "orange")

# GEA significance test
# showing the K = 5 estimated factors
plot(mod.lfmm2@U, pch = 19,
     xlab = "Latent Factors",
     ylab = "MeanPrecipitation During Warmest Quarter",
     col = as.factor(meta2$region))
legend("bottomleft", legend= levels(as.factor(meta2$region)), pch=16, col=unique(as.factor(meta2$region)))

# Simulate a matrix containing haploid genotypes
Y <- tcrossprod(as.matrix(X), B) +
  tcrossprod(U, V) +
  matrix(rnorm(n*L, sd = .5), nrow = n)
Y <- matrix(as.numeric(Y > 0), ncol = L)

# Fitting an LFMM with K = 5 factors
mod <- lfmm2(Y, PrecWmQ_env, K = 5)

# Computing P-values and plotting their minus log10 values
pv <- lfmm2.test(mod,
                 Y,
                 PrecWmQ_env,
                 linear = TRUE)
plot(-log10(pv$pvalues), col = "grey", cex = .6, pch = 19)
points(target, -log10(pv$pvalues[target]), col = "red")
