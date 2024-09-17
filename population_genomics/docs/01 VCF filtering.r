library(vcfR)

setwd("/gpfs1/cl/pbio3990/PopulationGenomics/")

list.files("variants/")

vcf <- read.vcfR("variants/Centaurea_filtered.vcf.gz")
vcf
head(vcf)

dna <- ape::read.dna("reference/GCA_030169165.1_ASM3016916v1_genomic.fa.gz", format="fasta")

gff <- read.table("reference/GCA_030169165.1_ASM3016916v1_genomic.fa.gz", sep="\t", quote="")

chr1 <- create.chromR(name="Chromosome 1", vcf=vcf, seq=dna, ann=gff)

plot(chr1)

#pdf(file="~/projects/eco_genomics/population_genomics/figures/ChromoPlot.pdf")
chromoqc(chr1, xlim=c(1e1, 1.1e8))
dev.off

DP <- extract.gt(vcf, element = "DP", as.numeric = T)

dim(DP)

DP[1:5,1:10]

quantile(DP)

#set zeros to NA
DP[DP==0] <- NA

#remove NA values
quantile(DP, na.rm = T)

# Visualize the matrix of depth and missingness in our VCF file

heatmap.bp(DP)
#subset 
heatmap.bp(DP[1:1000,],rlabels=F,clabels=F)

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
                                      cutoff=0.75)

vcf.filt.indMiss <- filter_biallelic(vcf.filt.indMiss)
vcf.filt.indMiss <- min_mac(vcf.filt.indMiss, min.mac = 1)

vcf.filt.indSNPMiss <- missing_by_snp(vcf.filt.indMiss, cutoff=0.5)

DP2 <- extract.gt(vcf.filt.indSNPMiss,
                  element="DP",
                  as.numeric = T)

heatmap.bp(DP2[5001:10000,],
           rlabels=F, clabels=F)

write.vcf(vcf.filt.indSNPMiss,
          "~/projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.vcf.gz")
