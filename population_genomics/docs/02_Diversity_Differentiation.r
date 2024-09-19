#Estimating diversity and genetic differentiation in the filtered Centaurea data

library(vcfR)
library(tidyverse)
library(qqman)

# helps solve plotting issues
X11.options(type="cairo")

#read in vcf file from our repo outputs/ directory

vcf<- read.vcfR("~/projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.vcf.gz")

#read in our metadata -- info on population of origin, what region pops come from, what continent, etc
meta <- read.csv("/gpfs1/cl/pbio3990/PopulationGenomics/metadata/meta4vcf.csv")

head(meta)#vcf has 593 samples
dim(meta)#meta has 629 samples

meta2 <- meta[meta$id %in% colnames(vcf@gt[,-1]),]

dim(meta2)

# calculate diversity stats using the genetic_diff function in vcfR
vcf.div <- genetic_diff(vcf,
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
          ylab="Fst among regions",
          suggestiveline = quantile(vcf.div.MHplot$Gst, 0.999))
