library(vcfR)
library(tidyverse)
library(qqman)
library(SNPfiltR)
library(LEA)

#new thresholds = 0.50 and 0.90

options(bitmapType = "cairo")

setwd("/gpfs1/cl/pbio3990/PopulationGenomics/")

vcf <- read.vcfR("variants/Centaurea_filtered.vcf.gz")

dna <- ape::read.dna("reference/GCA_030169165.1_ASM3016916v1_genomic.fa.gz", format="fasta")

gff <- read.table("reference/GCA_030169165.1_ASM3016916v1_genomic.fa.gz", sep="\t", quote="")

DP <- extract.gt(vcf, element = "DP", as.numeric = T)

DP[1:5,1:10]

quantile(DP)

#set zeros to NA
DP[DP==0] <- NA

#remove NA values
quantile(DP, na.rm = T)

# Visualize the matrix of depth and missingness in our VCF file

#heatmap.bp(DP)
#subset 
#heatmap.bp(DP[1:1000,],rlabels=F,clabels=F)

library(SNPfiltR)

#assign depth value and filter to a plot, addresses low depth problem here
vcf.filt <- hard_filter(vcf, depth=3) #could explore other depth values DP = 5, 10,...

#high depth problem
max_depth(vcf.filt)
#rule out values above twice the average as a general rule of thumb

vcf.filt <- max_depth(vcf.filt, maxdepth = 60)#filter out genotypes with >60 reads/SNP

#addressing missingness
meta <- read.csv("metadata/meta4vcf.csv", header=T)

meta2 <- meta[,c(1,4)]

names(meta2) <- c("id","pop")
meta2$id= as.factor(meta2$id)
meta2$pop=as.factor(meta2$pop)

#individual level missingness
vcf.filt.indMiss <- missing_by_sample(vcf.filt,
                                      popmap=meta2,
                                      cutoff=0.90)#change cutoff here!

vcf.filt.indMiss <- filter_biallelic(vcf.filt.indMiss)
vcf.filt.indMiss <- min_mac(vcf.filt.indMiss, min.mac = 1)

vcf.filt.indSNPMiss <- missing_by_snp(vcf.filt.indMiss, cutoff=0.90)#and here!

DP2 <- extract.gt(vcf.filt.indSNPMiss,
                  element="DP",
                  as.numeric = T)

#heatmap.bp(DP2[5001:10000,],
           #rlabels=F, clabels=F)

write.vcf(vcf.filt.indSNPMiss,
          "~/projects/eco_genomics/population_genomics/homework1/vcf_final.filtered_0.90missingness.vcf.gz")

#pt2
vcf<- read.vcfR("~/projects/eco_genomics/population_genomics/homework1/vcf_final.filtered_0.90missingness.vcf.gz")

meta <- read.csv("/gpfs1/cl/pbio3990/PopulationGenomics/metadata/meta4vcf.csv")

meta2 <- meta[meta$id %in% colnames(vcf@gt[,-1]),]

vcf.div <- genetic_diff(vcf,
                        pops=as.factor(meta2$region),
                        method = "nei")

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
          ylab="Fst among regions",
          suggestiveline = quantile(vcf.div.MHplot$Gst, 0.999))

write.csv(vcf.div.MHplot, "~/projects/eco_genomics/population_genomics/homework1/Genetic_Diff_byRegion_0.90missingness_hw1.csv",
          quote=F,
          row.names=F)

#make into long format and plot
vcf.div.MHplot %>% 
  as_tibble() %>%
  pivot_longer(c(4:9)) %>%
  ggplot(aes(x=value, fill=name)) +
  geom_histogram(position = "identity", alpha=0.5, bins=50) +
  labs(title="Genome-wide expected heterozygosity (Hs) 0.90", fill="Regions",
       x="Gene diversity within Regions", y= "Counts of SNPs")

#ggsave saves last plot
ggsave("Histogram_GenomeDiversity_byRegion_0.90missingness_hw1.pdf",
       path="~/projects/eco_genomics/population_genomics/homework1/")

#genome wide diversity, Hs, and StdDev
vcf.div.MHplot %>% 
  as_tibble() %>%
  pivot_longer(c(4:9)) %>%
  group_by(name) %>%
  filter(value!=0 & value<0.25) %>%
  summarise(avg_Hs = mean(value),StdDev_Hs = sd(value), N_Hs = n())

#pt3

# rename files for this homework

vcf.thin <- distance_thin(vcf, min.distance = 500)

meta <- read.csv("/gpfs1/cl/pbio3990/PopulationGenomics/metadata/meta4vcf.csv")

meta2 <- meta[meta$id %in% colnames(vcf@gt[, -1]),]

write.vcf(vcf.thin, "~/projects/eco_genomics/population_genomics/homework1/vcf_final.filtered.thinned_0.90missingness.vcf.gz")

#hide the uncompressed vcf file too big for github outside of our repo
system("gunzip -c ~/projects/eco_genomics/population_genomics/homework1/vcf_final.filtered.thinned_0.90missingness.vcf.gz > ~/vcf_final.filtered.thinned_0.90missingness.vcf")

geno <- vcf2geno(input.file = "~/vcf_final.filtered.thinned_0.90missingness.vcf",
                 output.file = "homework1/vcf_final.filtered.thinned_0.90missingness.geno")

CentPCA <- LEA::pca("homework1/vcf_final.filtered.thinned_0.90missingness.geno", scale=TRUE)

ggplot(as.data.frame(CentPCA$projections),
       aes(x=V1, y=V2, color=meta2$region, shape=meta2$continent))+
  geom_point(alpha=0.5) +
  labs(title="Centaurea genetic PCA",x="PC1", y="PC2", color="Region", shape="Continent")

ggsave("outputs/homework1/CentPCA_PC1vPC2_hw1.pdf", width=6, height=6, units="in")
