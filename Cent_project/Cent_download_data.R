# Required packages
# Loading worldclim/cimp6 bioclimatic data
library(terra)
library(geodata)
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


meta <- read.csv("/users/d/k/dkaupu/projects/eco_genomics/group_project/outputs/metafinal.csv", row.names = "X")

dim(meta)

meta2 <- meta[meta$id %in% colnames(vcf@gt[, -1]),]

dim(meta2)

#########
#Visualize coordinates
#########

genotype = LEA::read.geno("/users/l/e/lericsso/projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.thinned.geno")
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

geno_lfmm <- geno2lfmm(input.file = "/users/l/e/lericsso/projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.thinned.geno", output.file = "/users/l/e/lericsso/projects/eco_genomics/Cent_project/geno.lfmm", force = TRUE)

write.env(meta2[,9:11], output.file = "/users/l/e/lericsso/projects/eco_genomics/Cent_project/climate.env")
meta_env <- read.env("/users/l/e/lericsso/projects/eco_genomics/Cent_project/climate.env")

# main options
# K = number of ancestral populations
# entropy = TRUE computes the cross-entropy criterion,
# CPU = 4 is the number of CPU used (hidden input)
project = NULL
project = snmf("/users/l/e/lericsso/projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.thinned.geno",
               K = 1:10,
               entropy = TRUE,
               repetitions = 10,
               project = "new")

plot(project, col = "blue", pch = 19, cex = 1.2)

########
best = which.min(cross.entropy(project, K = 5))
my.colors <- c("tomato", "lightblue",
               "olivedrab", "gold")
barchart(project, K = 5, run = best,
         border = NA, space = 0,
         col = my.colors,
         xlab = "Individuals",
         ylab = "Ancestry proportions",
         main = "Ancestry matrix") -> bp
axis(1, at = 1:length(bp$order),
     labels = bp$order, las=1,
     cex.axis = .4)

######
p = snmf.pvalues(project,
                 entropy = TRUE,
                 ploidy = 2,
                 K = 4)
pvalues = p$pvalues
par(mfrow = c(2,1))
hist(pvalues, col = "orange")
plot(-log10(pvalues), pch = 19, col = "blue", cex = .5)

#######

project.missing = snmf("/users/l/e/lericsso/projects/eco_genomics/Cent_project/geno.lfmm", K = 4,
                       entropy = TRUE, repetitions = 10,
                       project = "new")

# select the run with the lowest cross-entropy value
best = which.min(cross.entropy(project.missing, K = 4))
# Impute the missing genotypes
impute(project.missing, "/users/l/e/lericsso/projects/eco_genomics/Cent_project/geno.lfmm",
       method = 'mode', K = 4, run = best)
## Missing genotype imputation for K = 4
## Missing genotype imputation for run = 4
## Results are written in the file: genoM.lfmm_imputed.lfmm
# Proportion of correct imputation results
dat.imp = read.lfmm("/users/l/e/lericsso/projects/eco_genomics/Cent_project/geno.lfmm_imputed.lfmm")
mean( tutorial.R[dat == 9] == dat.imp[dat == 9] )

######
geno_vcf <- read.vcfR("/users/l/e/lericsso/projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.thinned.vcf.gz")

mod.lfmm2 <- lfmm2(input = geno_vcf, env = meta_env, K = 2)
