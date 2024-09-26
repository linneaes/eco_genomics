library(tidyverse)
library(vcfR)
library(SNPfiltR)
library(LEA)

options(bitmapType = "cairo")

setwd("~/projects/eco_genomics/population_genomics/")

vcf <- read.vcfR("outputs/vcf_final.filtered.vcf.gz")

# We need to thin the SNPs for LD(linkage disequilibrium) before we run 
#PCA and Admixture analyses to satisfy the assumptions of independence among loci

vcf.thin <- distance_thin(vcf, min.distance = 500)

meta <- read.csv("/gpfs1/cl/pbio3990/PopulationGenomics/metadata/meta4vcf.csv")

dim(meta)

meta2 <- meta[meta$id %in% colnames(vcf.thin@gt[, -1]),]

dim(meta2)

write.vcf(vcf.thin, "outputs/vcf_final.filtered.thinned.vcf.gz")

#hide the uncompressed vcf file too big for github outside of our repo
system("gunzip -c ~/projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.thinned.vcf.gz > ~/vcf_final.filtered.thinned.vcf")

geno <- vcf2geno(input.file = "/gpfs1/home/l/e/lericsso/vcf_final.filtered.thinned.vcf",
                 output.file = "outputs/vcf_final.filtered.thinned.geno")

CentPCA <- LEA::pca("outputs/vcf_final.filtered.thinned.geno", scale=TRUE)

#You can load in the PCA results if you've done it previously
CentPCA <-load.pcaProject("vcf_final.filtered.thinned.pcaProject")

show(CentPCA)

plot(CentPCA)

#plot(CentPCA$projections, 
     #col=as.factor(meta2$region))
#legend("bottomright", legend=as.factor(unique(meta2$region)), 
                                       #fill=as.factor(unique(meta2$region)))

ggplot(as.data.frame(CentPCA$projections),
       aes(x=V1, y=V2, color=meta2$region, shape=meta2$continent))+
       geom_point(alpha=0.5) +
  labs(title="Centaurea genetic PCA",x="PC1", y="PC2", color="Region", shape="Continent")

ggsave("figures/CentPCA_PC1vPC2.pdf", width=6, height=6, units="in")

#PC2 v PC3
ggplot(as.data.frame(CentPCA$projections),
       aes(x=V2, y=V3, color=meta2$region, shape=meta2$continent))+
  geom_point(alpha=0.5) +
  labs(title="Centaurea genetic PCA",x="PC2", y="PC3", color="Region", shape="Continent")
