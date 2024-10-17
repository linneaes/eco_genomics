#Making contrasts between groups
library(eulerr)

#set up groups within DESeq object
dds$group <- factor(paste0(dds$DevTemp, dds$FinalTemp))
design(dds) <- ~ group
dds <- DESeq(dds)
dim(dds)
resultsNames(dds)
#[1] "Intercept"               "group_D18A33_vs_D18A28"  "group_D18BASE_vs_D18A28" "group_D22A28_vs_D18A28"  "group_D22A33_vs_D18A28" 
#[6] "group_D22BASE_vs_D18A28"

# 1. Compare baseline gene expression between developmental treatment groups
res_D18_BASE_D22_BASE <- results(dds, contrast=c("group", "D18BASE", "D22BASE"), alpha = 0.05)
res_D18_BASE_D22_BASE <- res_D18_BASE_D22_BASE[!is.na(res_D18_BASE_D22_BASE$padj),]
res_D18_BASE_D22_BASE <- res_D18_BASE_D22_BASE[order(res_D18_BASE_D22_BASE$padj),]
head(res_D18_BASE_D22_BASE)
summary(res_D18_BASE_D22_BASE)

# make a list of which genes are in our comparisons of interest are differentially expressed (list of DEGs)
degs_D18_BASE_D22_BASE <- row.names(res_D18_BASE_D22_BASE[res_D18_BASE_D22_BASE$padj <0.05, ])

plotMA(res_D18_BASE_D22_BASE, ylim=c(-4,4))

# 2. Compare gene expression between developmental temperature treatment groups at A28
res_D18_A28_D22_A28 <- results(dds, contrast=c("group", "D18A28", "D22A28"), alpha = 0.05)
res_D18_A28_D22_A28 <- res_D18_A28_D22_A28[!is.na(res_D18_A28_D22_A28$padj),]
res_D18_A28_D22_A28 <- res_D18_A28_D22_A28[order(res_D18_A28_D22_A28$padj),]
head(res_D18_A28_D22_A28)
summary(res_D18_A28_D22_A28)

# make a list of which genes are in our comparisons of interest are differentially expressed (list of DEGs)
degs_D18_A28_D22_A28 <- row.names(res_D18_A28_D22_A28[res_D18_A28_D22_A28$padj <0.05, ])

plotMA(res_D18_A28_D22_A28, ylim=c(-4,4))

# 3. Compare gene expression between developmental temperature treatment groups at A33
res_D18_A33_D22_A33 <- results(dds, contrast=c("group", "D18A33", "D22A33"), alpha = 0.05)
res_D18_A33_D22_A33 <- res_D18_A33_D22_A33[!is.na(res_D18_A33_D22_A33$padj),]
res_D18_A33_D22_A33 <- res_D18_A33_D22_A33[order(res_D18_A33_D22_A33$padj),]
head(res_D18_A33_D22_A33)
summary(res_D18_A33_D22_A33)

# make a list of which genes are in our comparisons of interest are differentially expressed (list of DEGs)
degs_D18_A33_D22_A33 <- row.names(res_D18_A33_D22_A33[res_D18_A33_D22_A33$padj <0.05, ])

plotMA(res_D18_A33_D22_A33, ylim=c(-4,4))



# Counts of differentially expressed genes between groups at each level
length(degs_D18_BASE_D22_BASE) #1935 differentially expressed genes between DevTemp groups at baseline
length(degs_D18_A28_D22_A28) #296
length(degs_D18_A33_D22_A33) #78

#look at overlaps in which genes are differentially expressed in multiple contrasts
length(intersect(degs_D18_BASE_D22_BASE, degs_D18_A28_D22_A28)) #107 degs overlapping between BASE and A28
length(intersect(degs_D18_BASE_D22_BASE, degs_D18_A33_D22_A33)) #44 degs overlapping btwn BASE and A33
length(intersect(degs_D18_A28_D22_A28, degs_D18_A33_D22_A33)) #29 degs overlapping btwen A28 and A33

length(intersect(degs_D18_BASE_D22_BASE, intersect(degs_D18_A28_D22_A28, degs_D18_A33_D22_A33))) #23 overlapping between all of them
