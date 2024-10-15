# Transcriptomics

## 10-8-2024 DESeq

We started a new unit on transcriptomics, beginning to compare the trascripts of Copepods grown at different temperatures when exposed to varying temperatures.

Developmental plasticity: same genes, different environment -\> different phenotypes

Experimental Questions: How do animals acclimate?: Physiological mechanisms of developmental plasticity: 1. Does the temp that they experience growing up effect UGT?, 2. How does gene expression response differ btwn temps? 3. What genes are expressed during development?

Factors: Development (levels = 18 degrees and 22 degrees C) and Final Temperature (levels = Control, 28, and 33 degrees C)

## 10-10-2024 DESeq2

We analyzed the Copepod data using DESeq and plotted the variation between treatments on a PCA plot.

Overdispersion = where variance \> mean -\> negative binomial distribution takes this into account

Gene Expression: explanatory variables = treatment conditions (DevTemp and FinalTemp)

colSums sums up all the columns, rowSums sums the rows

cex.names = size of the names

las = orientation of names

abline(h=mean(colSums(countsTableRound)), col="blue4", lwd=2) : adds a line on the barplot, h says where, lwd is the line width

nrow() number of rows

command - (shortcut to \<-)

command shift c (comments out a line)

## 10-15 Transcriptomics

We loaded in our day 1 script and data and then ordered it by significance. Then we looked at counts for a specific gene and compared the gene expression from the two developmental temperatures (18 and 22) in a plot.

We then made an MA plot to compare log fold change of differential gene expression vs. avg. gene expression. The very far right of the graph represents constitutely expressed/housekeeping genes that are always on. We see the expression of overdiversion, with a lot more variability than expected

We then made a volcano plot and saw more significance in upregulated genes in D22 (which is also the downregulated genes in D18). There had to be at least a doubling of expression to be considered significant.

We then made a heatmap to visualize the comparison of expression for different developmental and final temperatures. Each column is a sample, each row is a gene, and the color refers to magnitude of expression. We see some genes with high contrast in expression across samples, and a few genes that are significantly lowly expressed in the D18 individuals. The intensity of expression varies in the genes that all show up with the same color.

log2FoldChange = measure of difference in gene expression in units of doublings (instead of x10 in log10)
