#load in day 1 transcriptomics script to get conditions and counts table and run DESeq object

library(pheatmap)

options(bitmapType = "cairo")

resultsNames(dds)
##[1] "Intercept"             "DevTemp_D22_vs_D18"    "FinalTemp_A33_vs_A28"  "FinalTemp_BASE_vs_A28"

#pull out the results for Developmental Temperature 22 vs 18 
res_D22vsD18 <- results(dds, name = "DevTemp_D22_vs_D18", alpha =0.05)

#order by significance
res_D22vsD18 <- res_D22vsD18[order(res_D22vsD18$padj),]
head(res_D22vsD18)

#look at counts of a specific top gene that we're interested in to validate that the model is working
d <- plotCounts(dds, gene="TRINITY_DN140616_c0_g2_i1", int=(c("DevTemp", "FinalTemp")), returnData = TRUE)
d

p <- ggplot(d, aes(x=DevTemp, y=count, color=DevTemp, shape= FinalTemp))+
  theme_minimal() + theme(text=element_text(size=20), panel.grid.major = element_line(colour="grey"))
p <- p + geom_point(position = position_jitter(w=0.2,h=0), size=3)
p

# MA plot, log fold change of differential gene expression vs. avg. gene expression
plotMA(res_D22vsD18, ylim=c(-4,4))

# Volcano plot
# convert our DESeq results object into a data frame to plot
res_df <- as.data.frame(res_D22vsD18)

# add a column to dataframe to denote whether a gene is significantly differentially expressed or not
res_df$Significant <- ifelse(res_df$padj <0.05 & abs(res_df$log2FoldChange) > 1, "Significant", "Not Significant")

# plot
ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = Significant)) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = c("slateblue", "tomato")) +
  labs(x = "Log2 Fold Change", y = "log10 Adjusted P-value", title = "Volcano Plot")+
  theme_minimal()+
  theme(legend.position = "top") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "orange")+
  geom_vline(xintercept = c(-1, 1), linetype="dashed", color ="orange")
  
# Heat map
vsd <- vst(dds, blind = FALSE)

topgenes <- head(rownames(res_D22vsD18), 20)
mat <- assay(vsd)[topgenes, ]
df <- as.data.frame(colData(dds)[,c("DevTemp", "FinalTemp")])
pheatmap(mat, annotation_col=df, show_rownames=FALSE, cluster_cols=T, cluster_rows=T)

