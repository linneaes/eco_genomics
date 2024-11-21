# Required packages
# displaying images and maps
library(fields)
library(maps)
# Adjusting genotype-environment association models
library(LEA)
library(tidyverse)
library(vcfR)
library(SNPfiltR)
library(LEA)

options(bitmapType = "cairo")

setwd("~/projects/eco_genomics/population_genomics/")

vcf <- read.vcfR("~/projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.vcf.gz")

#system("gunzip -c ~/projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.vcf.gz > ~/vcf_final.filtered.vcf")

#geno <- vcf2geno(input.file = "/gpfs1/home/l/e/lericsso/vcf_final.filtered.vcf",
                 #output.file = "outputs/vcf_final.filtered.geno")

geno_vcf <- read.vcfR("/users/l/e/lericsso/projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.vcf.gz")

meta <- read.csv("/users/d/k/dkaupu/projects/eco_genomics/group_project/outputs/metafinal19.csv", row.names = "X")

dim(meta)

meta2 <- meta[meta$id %in% colnames(vcf@gt[, -1]),]

dim(meta2)

#########
#Visualize coordinates
#########

genotype = LEA::read.geno("/users/l/e/lericsso/projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.geno")
coordinates = as.matrix(read.table(meta2[,7:8]))
latitude <- meta2$latitude
longitude <- meta2$longitude

plot(longitude, latitude, cex = .4, col = "darkblue",
     xlab = "Longitude", ylab = "Latitude",
     main = "coordinates", las = 1)
maps::map(add = TRUE, interior = FALSE, col = "grey40")

#######
#climate data
#######

TempM <- meta2$TempM
TempR <- meta2$TempR
Prec <- meta2$Prec

geno_lfmm <- geno2lfmm(input.file = "/users/l/e/lericsso/projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.geno", output.file = "/users/l/e/lericsso/projects/eco_genomics/Cent_project/geno.lfmm", force = TRUE)

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
project = snmf("/users/l/e/lericsso/projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.geno",
               K = 1:10,
               entropy = TRUE,
               repetitions = 10,
               project = "new")

plot(project, col = "blue", pch = 19, cex = 1.2)

########
best = which.min(cross.entropy(project, K = 5))
library(RColorBrewer)
cols <- brewer.pal(6, "Set3")
cols
pal <- colorRampPalette(cols)
######
p = snmf.pvalues(project,
                 entropy = TRUE,
                 ploidy = 2,
                 K = 5)
pvalues = p$pvalues
par(mfrow = c(2,1))
hist(pvalues, col = "orange")
plot(-log10(pvalues), pch = 19, col = "blue", cex = .5)

#######

project.missing = snmf("/users/l/e/lericsso/projects/eco_genomics/Cent_project/geno.lfmm", K = 5,
                       entropy = TRUE, repetitions = 10,
                       project = "new")

# select the run with the lowest cross-entropy value
best = which.min(cross.entropy(project.missing, K = 5))
# Impute the missing genotypes
impute(project.missing, "/users/l/e/lericsso/projects/eco_genomics/Cent_project/geno.lfmm",
       method = 'mode', K = 5, run = best)
## Missing genotype imputation for K = 5
## Missing genotype imputation for run = 5
## Results are written in the file: /users/l/e/lericsso/projects/eco_genomics/Cent_project/geno.lfmm_imputed.lfmm
# Proportion of correct imputation results
dat.imp = read.lfmm("/users/l/e/lericsso/projects/eco_genomics/Cent_project/geno.lfmm_imputed.lfmm")
#mean( tutorial.R[dat == 9] == dat.imp[dat == 9] )

######
# LFMM2 Analysis Mean Temp
######

Y <- dat.imp
X <- meta2$TempM

mod.lfmm2 <- lfmm2(Y, X, K = 5)

# Simulate non-null effect sizes for 10 target loci
#individuals
n = 593
#loci
L = 3643

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
                 dat.imp,
                 TempM_env,
                 full = TRUE)
plot(-log10(pv$pvalues), col = "grey", cex = .5, pch = 19)
abline(h = -log10(0.1/510), lty = 2, col = "orange")

# GEA significance test
# showing the K = 5 estimated factors
plot(mod.lfmm2@U, pch = 19,
     xlab = "Latent Factors",
     ylab = "Mean Temperature",
     col = pal(as.factor(meta2$region)))
legend("bottomleft", legend= levels(as.factor(meta2$region)), pch=16, col=pal(as.factor(meta2$region)))

# Simulate a matrix containing haploid genotypes
Y <- tcrossprod(as.matrix(X), B) +
  tcrossprod(U, V) +
  matrix(rnorm(n*L, sd = .5), nrow = n)
Y <- matrix(as.numeric(Y > 0), ncol = L)

# Fitting an LFMM with K = 3 factors
mod <- lfmm2(Y, TempM_env, K = 5)

# Computing P-values and plotting their minus log10 values
pv <- lfmm2.test(mod,
                 Y,
                 TempM_env,
                 linear = TRUE)
plot(-log10(pv$pvalues), col = "grey", cex = .6, pch = 19)
points(target, -log10(pv$pvalues[target]), col = "red")

#########
# LFMM2 Temp range (Temp R)
#######
Y <- dat.imp
X <- meta2$TempR

mod.lfmm2 <- lfmm2(Y, X, K = 5)

# Simulate non-null effect sizes for 10 target loci
#individuals
n = 593
#loci
L = 3643
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
                 dat.imp,
                 TempR_env,
                 full = TRUE)
plot(-log10(pv$pvalues), col = "grey", cex = .5, pch = 19)
abline(h = -log10(0.1/510), lty = 2, col = "orange")

# GEA significance test
# showing the K = 5 estimated factors
plot(mod.lfmm2@U, pch = 19,
     xlab = "Latent Factors",
     ylab = "Temperature Range",
     col = pal(as.factor(meta2$region)))
legend("bottomleft", legend= levels(as.factor(meta2$region)), pch=16, col=pal(as.factor(meta2$region)))

# Simulate a matrix containing haploid genotypes
Y <- tcrossprod(as.matrix(X), B) +
  tcrossprod(U, V) +
  matrix(rnorm(n*L, sd = .5), nrow = n)
Y <- matrix(as.numeric(Y > 0), ncol = L)

# Fitting an LFMM with K = 3 factors
mod <- lfmm2(Y, TempR_env, K = 5)

# Computing P-values and plotting their minus log10 values
pv <- lfmm2.test(mod,
                 Y,
                 TempR_env,
                 linear = TRUE)
plot(-log10(pv$pvalues), col = "grey", cex = .6, pch = 19)
points(target, -log10(pv$pvalues[target]), col = "red")


#######
# LFMM2 Prec
#######
Y <- dat.imp
X <- meta2$Prec

mod.lfmm2 <- lfmm2(Y, X, K = 5)

# Simulate non-null effect sizes for 10 target loci
#individuals
n = 593
#loci
L = 3643
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
                 dat.imp,
                 Prec_env,
                 full = TRUE)
plot(-log10(pv$pvalues), col = "grey", cex = .5, pch = 19)
abline(h = -log10(0.1/510), lty = 2, col = "orange")

# GEA significance test
# showing the K = 5 estimated factors
plot(mod.lfmm2@U, pch = 19,
     xlab = "Latent Factors",
     ylab = "Precipitation",
     col = pal(as.factor(meta2$region)))
legend("bottomleft", legend= levels(as.factor(meta2$region)), pch=16, col=pal(as.factor(meta2$region)))

# Simulate a matrix containing haploid genotypes
Y <- tcrossprod(as.matrix(X), B) +
  tcrossprod(U, V) +
  matrix(rnorm(n*L, sd = .5), nrow = n)
Y <- matrix(as.numeric(Y > 0), ncol = L)

# Fitting an LFMM with K = 3 factors
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
X <- meta2$MDiurnalR

mod.lfmm2 <- lfmm2(Y, X, K = 5)

# Simulate non-null effect sizes for 10 target loci
#individuals
n = 593
#loci
L = 3643

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
                 dat.imp,
                 MDiurnalR_env,
                 full = TRUE)
plot(-log10(pv$pvalues), col = "grey", cex = .5, pch = 19)
abline(h = -log10(0.1/510), lty = 2, col = "orange")

# GEA significance test
# showing the K = 5 estimated factors
plot(mod.lfmm2@U, pch = 19,
     xlab = "Latent Factors",
     ylab = "Mean Diurnal Range",
     col = pal(as.factor(meta2$region)))
legend("bottomleft", legend= levels(as.factor(meta2$region)), pch=16, col=pal(as.factor(meta2$region)))

# Simulate a matrix containing haploid genotypes
Y <- tcrossprod(as.matrix(X), B) +
  tcrossprod(U, V) +
  matrix(rnorm(n*L, sd = .5), nrow = n)
Y <- matrix(as.numeric(Y > 0), ncol = L)

# Fitting an LFMM with K = 3 factors
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
X <- meta2$PrecDM

mod.lfmm2 <- lfmm2(Y, X, K = 5)

# Simulate non-null effect sizes for 10 target loci
#individuals
n = 593
#loci
L = 3643

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
                 dat.imp,
                 PrecDM_env,
                 full = TRUE)
plot(-log10(pv$pvalues), col = "grey", cex = .5, pch = 19)
abline(h = -log10(0.1/510), lty = 2, col = "orange")

# GEA significance test
# showing the K = 5 estimated factors
plot(mod.lfmm2@U, pch = 19,
     xlab = "Latent Factors",
     ylab = "Precipitation Driest Month",
     col = pal(as.factor(meta2$region)))
legend("bottomleft", legend= levels(as.factor(meta2$region)), pch=16, col=pal(as.factor(meta2$region)))

# Simulate a matrix containing haploid genotypes
Y <- tcrossprod(as.matrix(X), B) +
  tcrossprod(U, V) +
  matrix(rnorm(n*L, sd = .5), nrow = n)
Y <- matrix(as.numeric(Y > 0), ncol = L)

# Fitting an LFMM with K = 3 factors
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
X <- meta2$PrecWmQ

mod.lfmm2 <- lfmm2(Y, X, K = 5)

# Simulate non-null effect sizes for 10 target loci
#individuals
n = 593
#loci
L = 3643

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
                 dat.imp,
                 PrecWmQ_env,
                 full = TRUE)
plot(-log10(pv$pvalues), col = "grey", cex = .5, pch = 19)
abline(h = -log10(0.1/510), lty = 2, col = "orange")

# GEA significance test
# showing the K = 5 estimated factors
plot(mod.lfmm2@U, pch = 19,
     xlab = "Latent Factors",
     ylab = "Precipitation Warmest Quarter",
     col = pal(as.factor(meta2$region)))
legend("bottomleft", legend= levels(as.factor(meta2$region)), pch=16, col=pal(as.factor(meta2$region)))

# Simulate a matrix containing haploid genotypes
Y <- tcrossprod(as.matrix(X), B) +
  tcrossprod(U, V) +
  matrix(rnorm(n*L, sd = .5), nrow = n)
Y <- matrix(as.numeric(Y > 0), ncol = L)

# Fitting an LFMM with K = 3 factors
mod <- lfmm2(Y, PrecWmQ_env, K = 5)

# Computing P-values and plotting their minus log10 values
pv <- lfmm2.test(mod,
                 Y,
                 PrecWmQ_env,
                 linear = TRUE)
plot(-log10(pv$pvalues), col = "grey", cex = .6, pch = 19)
points(target, -log10(pv$pvalues[target]), col = "red")
