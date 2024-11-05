### Code for Question 3 Analysis

# load libraries
library(DESeq2)
library(ggplot2)
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


#################################################################################
# Contrasts
#################################################################################
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

length(intersect(degs_D18_BASE_D18_A28, degs_D18_BASE_D18_A33)) #34 degs overlapping btwen BASEs

41-34 = 7 #A28
332-34 = 298 #A33


myEuler <- euler(c("A28" = 7, "A33" = 298,
                   "A28&A33" = 34))

plot(myEuler, lty=1:3, quantities=TRUE, fill = c("goldenrod", "cornflowerblue", "olivedrab"))


#################################################################################
# DevTemp 22

# 1. Compare DevTemp 22 at Baseline and Acute 28
res_D22_BASE_D22_A28 <- results(dds, contrast=c("group", "D22BASE", "D22A28"), alpha = 0.05)
res_D22_BASE_D22_A28 <- res_D22_BASE_D22_A28[!is.na(res_D22_BASE_D22_A28$padj),]
res_D22_BASE_D22_A28 <- res_D22_BASE_D22_A28[order(res_D22_BASE_D22_A28$padj),]
head(res_D22_BASE_D22_A28)
summary(res_D22_BASE_D22_A28)

# make a list of which genes are in our comparisons of interest are differentially expressed (list of DEGs)
degs_D22_BASE_D22_A28 <- row.names(res_D22_BASE_D22_A28[res_D22_BASE_D22_A28$padj <0.05, ])

plotMA(res_D22_BASE_D22_A28, ylim=c(-4,4))



# 2. Compare DevTemp 22 at Baseline and Acute 33
res_D22_BASE_D22_A33 <- results(dds, contrast=c("group", "D22BASE", "D22A33"), alpha = 0.05)
res_D22_BASE_D22_A33 <- res_D22_BASE_D22_A33[!is.na(res_D22_BASE_D22_A33$padj),]
res_D22_BASE_D22_A33 <- res_D22_BASE_D22_A33[order(res_D22_BASE_D22_A33$padj),]
head(res_D22_BASE_D22_A33)
summary(res_D22_BASE_D22_A33)

# make a list of which genes are in our comparisons of interest are differentially expressed (list of DEGs)
degs_D22_BASE_D22_A33 <- row.names(res_D22_BASE_D22_A33[res_D22_BASE_D22_A33$padj <0.05, ])

plotMA(res_D22_BASE_D22_A33, ylim=c(-4,4))


# 3. Compare DevTemp 22 at Acute A28 and Acute 33
res_D22_A28_D22_A33 <- results(dds, contrast=c("group", "D22A28", "D22A33"), alpha = 0.05)
res_D22_A28_D22_A33 <- res_D22_A28_D22_A33[!is.na(res_D22_A28_D22_A33$padj),]
res_D22_A28_D22_A33 <- res_D22_A28_D22_A33[order(res_D22_A28_D22_A33$padj),]
head(res_D22_A28_D22_A33)
summary(res_D22_A28_D22_A33)

# make a list of which genes are in our comparisons of interest are differentially expressed (list of DEGs)
degs_D22_A28_D22_A33 <- row.names(res_D22_A28_D22_A33[res_D22_A28_D22_A33$padj <0.05, ])

plotMA(res_D22_A28_D22_A33, ylim=c(-4,4))


# Counts of differentially expressed genes between groups at each level
length(degs_D22_BASE_D22_A28) #289 diff. exp. genes between D22 at Baseline and A28
length(degs_D22_BASE_D22_A33) #1564 diff. exp. genes btwn D22 at Baseline and A33

length(intersect(degs_D22_BASE_D22_A28, degs_D22_BASE_D22_A33)) #144 degs overlapping btwen BASEs

289-144 = 145 #A28
1564-144 = 1420 #A33


myEuler <- euler(c("A28" = 145, "A33" = 1420,
                   "A28&A33" = 144))

plot(myEuler, lty=1:3, quantities=TRUE, fill = c("tomato", "slateblue2", "orchid"))

##################################################################################
# Scatter Plots
#################################################################################

# contrast D18_A28vsBASE
res_D18_BASEvsA28 <- as.data.frame(results(dds, contrast= c("group", "D18BASE", "D18A28"), alpha=0.05))

# contrast D18_A33vsBASE
res_D18_BASEvsA33 <- as.data.frame(results(dds, contrast= c("group", "D18BASE", "D18A33"), alpha=0.05))

# merge dataframes
res_df18 <- merge(res_D18_BASEvsA28, res_D18_BASEvsA33, by = "row.names", suffixes = c(".28", ".33"))

rownames(res_df18) <- res_df18$Row.names
res_df18 <- res_df18[ , -1] #get rid of first column

res_df18 <- res_df18 %>%
  mutate(fill = case_when(
    padj.28 < 0.05 & stat.28 < 0 ~ "goldenrod2",
    padj.28 < 0.05 & stat.28 > 0 ~ "cornflowerblue", 
    padj.33 < 0.05 & stat.33 < 0 ~ "tomato2", 
    padj.33 < 0.05 & stat.33 > 0 ~ "orchid"
  ))

# Count the number of points per fill color
color_counts <- res_df18 %>%
  group_by(fill) %>%
  summarise(count = n())

label_positions <- data.frame(
  fill = c("goldenrod2", "cornflowerblue", "tomato2", "orchid"),
  x_position = c(-3, 5, 1, 0),
  y_position = c(-8, 0, -5, 9)
)

label_data <- merge(color_counts, label_positions, by = "fill")

# Plot
plot18 <- ggplot(res_df18, aes(x= log2FoldChange.28, y= log2FoldChange.33, color = fill)) +
  geom_point(alpha = 0.8) +
  scale_color_identity() +
  geom_text(data = label_data, aes(x = x_position, y = y_position, label = count, color = fill),
            size = 5) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  geom_abline(intercept = 0, slope = -1, linetype = "dashed", color = "grey") +
  xlim(-10, 10) + ylim(-10, 10) +
  labs(x = "Log2FoldChange 28 vs BASE at 18", 
       y = "Log2FoldChange 33 vs BASE at 18",
       title = "How does response to 28 C and 33 C vary at DevTemp 18?") +
  theme_minimal()

plot18




# contrast D22_A28vsBASE
res_D22_BASEvsA28 <- as.data.frame(results(dds, contrast= c("group", "D22BASE", "D22A28"), alpha=0.05))

# contrast D22_A33vsBASE
res_D22_BASEvsA33 <- as.data.frame(results(dds, contrast= c("group", "D22BASE", "D22A33"), alpha=0.05))

# merge dataframes
res_df22 <- merge(res_D22_BASEvsA28, res_D22_BASEvsA33, by = "row.names", suffixes = c(".28", ".33"))

rownames(res_df22) <- res_df22$Row.names
res_df22 <- res_df22[ , -1] #get rid of first column

res_df22 <- res_df22 %>%
  mutate(fill = case_when(
    padj.28 < 0.05 & stat.28 < 0 ~ "goldenrod2",
    padj.28 < 0.05 & stat.28 > 0 ~ "cornflowerblue", 
    padj.33 < 0.05 & stat.33 < 0 ~ "tomato2", 
    padj.33 < 0.05 & stat.33 > 0 ~ "orchid"
  ))

# Count the number of points per fill color
color_counts <- res_df22 %>%
  group_by(fill) %>%
  summarise(count = n())

label_positions <- data.frame(
  fill = c("goldenrod2", "cornflowerblue", "tomato2", "orchid"),
  x_position = c(-2, 5, 1, 0),
  y_position = c(-8, 0, -5, 5)
)

label_data <- merge(color_counts, label_positions, by = "fill")

# Plot
plot22 <- ggplot(res_df22, aes(x= log2FoldChange.28, y= log2FoldChange.33, color = fill)) +
  geom_point(alpha = 0.8) +
  scale_color_identity() +
  geom_text(data = label_data, aes(x = x_position, y = y_position, label = count, color = fill),
            size = 5) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  geom_abline(intercept = 0, slope = -1, linetype = "dashed", color = "grey") +
  xlim(-10, 10) + ylim(-10, 10) +
  labs(x = "Log2FoldChange 28 vs BASE at 22", 
       y = "Log2FoldChange 33 vs BASE at 22",
       title = "How does response to 28 C and 33 C vary at DevTemp 22?") +
  theme_minimal()

plot22


combined_plot <- grid.arrange(plot18, plot22, ncol = 2)

ggsave("~/projects/eco_genomics/transcriptomics/Homework2/combined_scatter_plot_HW2.png", combined_plot, width = 12, height = 6)
