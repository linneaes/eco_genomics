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

meta2 <- meta[meta$id %in% colnames(vcf@gt[, -1]),]

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

# Now we will run admixture analyses and create plots
# For Admixture, we're going to use the LEA R package
# The function inside LEA is called "snmf"

CentAdmix <- snmf("outputs/vcf_final.filtered.thinned.geno", 
                  K=1:10,
                  entropy=T,
                  repetitions=3,
                  project="new") #if you're adding to this analysis later, you could choose project="continue"

par(mfrow= c(2,1))
plot(CentAdmix, col="blue4", main="SNMF") #this plots the cross-entropy score we can use for selecting models with K values that fit our data well
plot(CentPCA$eigenvalues[1:10], ylab="Eigenvalues", xlab= "number of PCs", col="blue4", main="PCA")
dev.off()

myK=5

CE = cross.entropy(CentAdmix, K=myK)
best = which.min(CE)

myKQ = Q(CentAdmix, K=myK, run=best)

#stitch together with metadata
myKQmeta = cbind(myKQ, meta2)

my.colors = c("blue4", "goldenrod","tomato", "lightblue", "olivedrab")

myKQmeta = as_tibble(myKQmeta) %>%
  group_by(continent) %>%
  arrange(region, pop, .by_group = TRUE)

pdf("figures/Admixture_K4.pdf", width=10, height=5)
barplot(as.matrix(t(myKQmeta[ , 1:myK])), 
        border=NA,
        space=0,
        col=my.colors[1:myK],
        xlab= "Geographic regions",
        ylab= "Ancestry proportions",
        main= paste0("Ancestry matrix K=",myK))
axis(1,
     at=1:length(myKQmeta$region),
     labels=myKQmeta$region,
     tick=F,
     cex.axis=0.5,
     las=3)
dev.off()
