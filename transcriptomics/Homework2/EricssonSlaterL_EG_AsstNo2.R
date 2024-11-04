### Code for Question 3 Analysis

# load libraries
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(eulerr)
library(dplyr)
library(tidyr)
library(gridExtra)

options(bitmapType="cairo")

setwd("~/projects/eco_genomics/transcriptomics/")

################################################################################
# DESeq2
# Import counts matrix

countsTable <- read.table("/gpfs1/cl/pbio3990/Transcriptomics/tonsa_counts.txt", 
                          header = TRUE, row.names = 1)

countsTableRound <- round(countsTable) # because DEseq2 doesn't like decimals
tail(countsTableRound) #see bottom of data

conds <- read.delim("/gpfs1/cl/pbio3990/Transcriptomics/experimental_details.txt",
                    header = TRUE, stringsAsFactors = TRUE, row.names = 1)
head(conds)

################################################################################
#
#Explore counts matrix
#
################################################################################

#let's see how many reads we have from each sample
colSums(countsTableRound)
mean(colSums(countsTableRound))

barplot(colSums(countsTableRound), names.arg=colnames(countsTableRound),
        cex.names = 0.5, las = 2, ylim = c(0,30000000))
abline(h=mean(colSums(countsTableRound)), col="blue4", lwd=2)

#the average number of counts per gene
rowSums(countsTableRound)
mean(rowSums(countsTableRound)) # mean = 3244.739
median(rowSums(countsTableRound)) # median = 64

#look at means of columns (2), and rows (1)
#gives a sense of variation in sequencing effort across samples
apply(countsTableRound, 2 , mean)
apply(countsTableRound, 1, mean)

################################################################################
#
#Start analysis in DESeq2
#
################################################################################


dds <- DESeqDataSetFromMatrix(countData= countsTableRound, colData = conds, 
                              design = ~ DevTemp + FinalTemp)
dim(dds)

dds<- dds[rowSums(counts(dds) >= 10) >= 15,] 
nrow(dds) # 35527 = number of transcripts with more than 10 reads and more than or equal to 15 samples

# run the DESeq model to test for global differential gene expression
dds<- DESeq(dds)

#list the results you've generated
resultsNames(dds)
# [1] "Intercept"             "DevTemp_D22_vs_D22"    "FinalTemp_A33_vs_A28"  "FinalTemp_BASE_vs_A28"

# visualize our global gene expression patterns using PCA
# first we need to transform the data for plotting using variance stabilization

vsd <- vst(dds, blind= FALSE)

pcaData <- plotPCA(vsd, intgroup= c("DevTemp", "FinalTemp"), returnData = TRUE)
percentVar <-  round(100*attr(pcaData, "percentVar"))

final_temp_colors <- c("BASE" = "grey", "A28" = "hotpink", "A33" = "red")
shapes_choose <- c("D22" = 16, "D22" = 22)

p <- ggplot(pcaData, aes(PC1, PC2, color=FinalTemp, shape= DevTemp))+
  geom_point(size=5)+
  scale_shape_manual(values = shapes_choose )+
  scale_color_manual(values = final_temp_colors) +
  labs(x = paste0('PC1: ', percentVar[1], ' %'),
       y= paste0('PC2: ', percentVar[2], ' %')) +
  theme_bw(base_size = 16)

p 

#################################################################################
# # Transcriptomics cont.: MA plot, Volcano plot, Heatmap
# 
# #pull out the results for Developmental Temperature 22 vs 22 
# res_D22vsD18 <- results(dds, name = "DevTemp_D22_vs_D18", alpha =0.05)
# 
# #order by significance
# res_D22vsD18 <- res_D22vsD18[order(res_D22vsD18$padj),]
# head(res_D22vsD18)
# 
# #look at counts of a specific top gene that we're interested in to validate that the model is working
# d <- plotCounts(dds, gene="TRINITY_DN140616_c0_g2_i1", int=(c("DevTemp", "FinalTemp")), returnData = TRUE)
# d
# 
# p <- ggplot(d, aes(x=DevTemp, y=count, color=DevTemp, shape= FinalTemp))+
#   theme_minimal() + theme(text=element_text(size=20), panel.grid.major = element_line(colour="grey"))
# p <- p + geom_point(position = position_jitter(w=0.2,h=0), size=3)
# p
# 
# # MA plot, log fold change of differential gene expression vs. avg. gene expression
# plotMA(res_D22vsD18, ylim=c(-4,4))
# 
# # Volcano plot
# # convert our DESeq results object into a data frame to plot
# res_df <- as.data.frame(res_D22vsD18)
# 
# # add a column to dataframe to denote whether a gene is significantly differentially expressed or not
# res_df$Significant <- ifelse(res_df$padj <0.05 & abs(res_df$log2FoldChange) > 1, "Significant", "Not Significant")
# 
# # plot
# ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = Significant)) +
#   geom_point(alpha = 0.8) +
#   scale_color_manual(values = c("slateblue", "tomato")) +
#   labs(x = "Log2 Fold Change", y = "log10 Adjusted P-value", title = "Volcano Plot")+
#   theme_minimal()+
#   theme(legend.position = "top") +
#   geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "orange")+
#   geom_vline(xintercept = c(-1, 1), linetype="dashed", color ="orange")
# 
# # Heat map
# vsd <- vst(dds, blind = FALSE)
# 
# topgenes <- head(rownames(res_D22vsD18), 20)
# mat <- assay(vsd)[topgenes, ]
# df <- as.data.frame(colData(dds)[,c("DevTemp", "FinalTemp")])
# pheatmap(mat, annotation_col=df, show_rownames=FALSE, cluster_cols=T, cluster_rows=T)

#################################################################################
# Contrasts

#set up groups within DESeq object
dds$group <- factor(paste0(dds$DevTemp, dds$FinalTemp))
design(dds) <- ~ group
dds <- DESeq(dds)
dim(dds)
resultsNames(dds)
#[1] "Intercept"               "group_D18A33_vs_D18A28"  "group_D18BASE_vs_D18A28"
#[4] "group_D22A28_vs_D18A28"  "group_D22A33_vs_D18A28"  "group_D22BASE_vs_D18A28"

#################################################################################
# DevTemp 18

# 1. Compare DevTemp 18 at Baseline and Acute 28
res_D18_BASE_D18_A28 <- results(dds, contrast=c("group", "D18BASE", "D18A28"), alpha = 0.05)
res_D18_BASE_D18_A28 <- res_D18_BASE_D18_A28[!is.na(res_D18_BASE_D18_A28$padj),]
res_D18_BASE_D18_A28 <- res_D18_BASE_D18_A28[order(res_D18_BASE_D18_A28$padj),]
head(res_D18_BASE_D18_A28)
summary(res_D18_BASE_D18_A28)

# make a list of which genes are in our comparisons of interest are differentially expressed (list of DEGs)
degs_D18_BASE_D18_A28 <- row.names(res_D18_BASE_D18_A28[res_D18_BASE_D18_A28$padj <0.05, ])

plotMA(res_D18_BASE_D18_A28, ylim=c(-4,4))



# 2. Compare DevTemp 18 at Baseline and Acute 33
res_D18_BASE_D18_A33 <- results(dds, contrast=c("group", "D18BASE", "D18A33"), alpha = 0.05)
res_D18_BASE_D18_A33 <- res_D18_BASE_D18_A33[!is.na(res_D18_BASE_D18_A33$padj),]
res_D18_BASE_D18_A33 <- res_D18_BASE_D18_A33[order(res_D18_BASE_D18_A33$padj),]
head(res_D18_BASE_D18_A33)
summary(res_D18_BASE_D18_A33)

# make a list of which genes are in our comparisons of interest are differentially expressed (list of DEGs)
degs_D18_BASE_D18_A33 <- row.names(res_D18_BASE_D18_A33[res_D18_BASE_D18_A33$padj <0.05, ])

plotMA(res_D18_BASE_D18_A33, ylim=c(-4,4))


# 3. Compare DevTemp 18 at Acute A28 and Acute 33
res_D18_A28_D18_A33 <- results(dds, contrast=c("group", "D18A28", "D18A33"), alpha = 0.05)
res_D18_A28_D18_A33 <- res_D18_A28_D18_A33[!is.na(res_D18_A28_D18_A33$padj),]
res_D18_A28_D18_A33 <- res_D18_A28_D18_A33[order(res_D18_A28_D18_A33$padj),]
head(res_D18_A28_D18_A33)
summary(res_D18_A28_D18_A33)

# make a list of which genes are in our comparisons of interest are differentially expressed (list of DEGs)
degs_D18_A28_D18_A33 <- row.names(res_D18_A28_D18_A33[res_D18_A28_D18_A33$padj <0.05, ])

plotMA(res_D18_A28_D18_A33, ylim=c(-4,4))


# Counts of differentially expressed genes between groups at each level
length(degs_D18_BASE_D18_A28) #41 diff. exp. genes between D18 at Baseline and A28
length(degs_D18_BASE_D18_A33) #332 diff. exp. genes btwn D18 at Baseline and A33
length(degs_D18_A28_D18_A33) #234 diff. exp. genes btwn D18 at A28 and A33

#look at overlaps in which genes are differentially expressed in multiple contrasts
length(intersect(degs_D18_BASE_D18_A28, degs_D18_A28_D18_A33)) #18 degs overlapping btwn BASE and A28
length(intersect(degs_D18_BASE_D18_A33, degs_D18_A28_D18_A33)) #163 degs overlapping btwn BASE and A33
length(intersect(degs_D18_BASE_D18_A28, degs_D18_BASE_D18_A33)) #34 degs overlapping btwen BASEs

length(intersect(degs_D18_A28_D18_A33, intersect(degs_D18_BASE_D18_A28, degs_D18_BASE_D18_A33))) #17 overlapping between all of them

# 
# # Calculate the number of unique genes in each portion of the Euler plot
# 1935-107-34+17 # 2207 genes differentially expressed uniquely between D18 at A

# 296-107-29+23 # 223 genes uniquely expressed when exposed to 28
# 78-44-29+23 # 28 genes uniquely expressed when exposed to 33
# 
# 107-23 # 84 unique to BASE and A28
# 44-23 # 21 unique to BASE and A33
# 29-23 # 6 unique to A28 and A33
# 
myEuler <- euler(c("BASE" = 2207, "A28" = 223, "A33" = 28,
                   "BASE&A28" = 84, "BASE&A33" = 21, "A28&A33" = 6,
                   "BASE&A28&A33" = 23))

plot(myEuler, lty=1:3, quantities=TRUE)


#####################################################################################
# DevTemp 22


