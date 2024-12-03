#Estimating diversity and genetic differentiation in the filtered Centaurea data

library(vcfR)
library(tidyverse)
library(qqman)

#helps with plotting errors on RStudio on VACC
options(bitmapType = "cairo")

vcf <- read.vcfR("~/projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.vcf.gz")

geno_vcf <- read.vcfR("/users/l/e/lericsso/projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.vcf.gz")

meta <- read.csv("/users/d/k/dkaupu/projects/eco_genomics/group_project/outputs/metafinal19.csv", row.names = "X")

dim(meta)

meta2 <- meta[meta$id %in% colnames(vcf@gt[, -1]),]

dim(meta2)

#############

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


############

# calculate diversity stats using the genetic_diff function in vcfR
vcf.div <- genetic_diff(geno_vcf,
                        pops=as.factor(meta2$region),
                        method = "nei")

str(vcf.div)

chr.main <- unique(vcf.div$CHROM)[1:8]

chrnum <- as.data.frame(cbind(chr.main, seq(1,8,1)))

vcf.div.MHplot <- left_join(chrnum,vcf.div, join_by(chr.main==CHROM))

vcf.div.MHplot <- vcf.div.MHplot %>%
  filter(Gst>0)  %>%  
  mutate(SNP=paste0(chr.main, "_", POS))  

vcf.div.MHplot$V2 = as.numeric(vcf.div.MHplot$V2)
vcf.div.MHplot$POS = as.numeric(vcf.div.MHplot$POS)


manhattan(vcf.div.MHplot,
          chr="V2",
          bp="POS",
          p="Gst",
          col=c("blue4","orange3"),
          logp=F,
          ylab="-log10pv$pvalues",
          suggestiveline = quantile(vcf.div.MHplot$Gst, 0.999))


