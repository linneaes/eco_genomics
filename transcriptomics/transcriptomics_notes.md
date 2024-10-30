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

## 10-17 Differentially Expressed Genes

Factors= DevTemp (D18 and D22) and Final Temp (Base, A28, and A33)

We combined the treatments into one variable called group, which had separate levels based on devTemp and FinalTemp. We're comparing gene expression between groups and within groups.

We started by looking at the counts of differentially expressed genes in each group at the baseline, then at A28 and A33. We did not find a large difference in the number of genes that were upregulated or downregulated.

We see a much larger number of differentially expressed genes at baseline, then less when exposed to 28, and even less at 33. The higher the acute exposure, the lower the difference between gene expression

We looked at overlap between differential gene expression at different acute temps vs the baseline. We found 23 differentially expressed genes overlapping between all 3 temps

ctrl F and select a region to replace words

## 10-22 Euler plots, Scatter plots, and WGCNA

We finished creating euler plots for the data to compare the overlap of differentially expressed genes between treatment groups (BASE, A28, A33).

We then created scatter plots of responses from different devTemps when exposed to A28 and A33. We contrasted the baseline treatment and A28 for the devtemp 18 group, and then merged the two dataframes. We saw a huge upregulation when exposed to 28 with the baseline at 22, and a downregulation from 18.

## 10-24

We finished the scatter plot for responses to A33 exposure and then placed the two plots together in a figure. This lets us compare all the differentially expressed genes at A28 vs A33. Many more genes were differentially expressed in A33, and there was more downregulation of the DevTemp 18 genes at A33. We also saw a large upregulation in genes at 22.

We then looked at WGCNA:

We looked at the good genes on a dendrogram and in a PCA plot to see the outliers in the 7 genes. There is one very strong outlier identified with both methods, but they are left in because they are within the expected variation. We then moved on to normalization, and then to network construction, where we looked at the scale free topology of the genes and their connectivity to other genes. We can use the plots that we made to choose what power will work well for our data (we chose about 26 because it just clears the 0.8 threshold and is slightly above the lowest connectivity value).

## 10-29-24 WGCNA continued

Power is how well topology of the the network represents something biological. "Scale free topology" - in a network you can have regions that are super connected with each other and some that are less tightly connected. Power = weakness of the correlation between genes -\> we used the plots to see which power represented the scale free top. model well and didn't let connectivity get too low
