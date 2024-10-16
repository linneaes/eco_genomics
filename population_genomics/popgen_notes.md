# Coding and data notes for the population genomics module

## Author: Linnea Ericsson Slater

### 09-10-2024 - Intro to Centaurea GBS data and working with VCF files

We'll be analyzing the GBS data from 3 regions (EU, NE, PNW) starting today with variant call format files (VCFs)

We learned how to set up our notes file and use markdown shortcuts. We then used the command line by opening VACC Shell Access and learned how to switch our directories over to the example data, then how to open it with zcat and use \| head to look at the first 10 lines. We then learned how to change the amount of data lines in the head by adding -n and a number (example: zcat file.gz \| head -n 4).

cd = change directory

zcat is used to view a zip file

\|s -\| = long list -\> can be abbreviated as ll

head = view first 10 lines of script

\| "pipe" = send output of one function to another

If you are in the correct directory w/ a specific file or directory, you can use tab to auto-complete the name

Our example data was a Q-file, so we saw the Q-scores in I's and G's according to what the score was, along with the sequence itself, and the barcode identifying the individual sample.

### 09-12-2024 - Viewing VCF files and talking about filtering

We learned how to open samtools on the command line and use it to look at a bam file to see the summary of reads and the sequencing depth from the data.

cd .. gets you out of a directory one time back

q to quit

### 09-17- VCF Filtering and diversity stats

We filtered our dataset to account for the 3 problems: Depth, Missingness, and Low-Frequency alleles. We also viewed the heatmaps of the data as we were trimming it down, which gave us a few figures to approximate our data and see what we might need to account for.

list.files("variants/")

```         
GT:PL:DP:AD
GT: genotype, 0/0:homozygote, 0/1:heterozygote, ./.:no data
PL: phred-scaled genotype likelihood
DP:Depth
AD:Allele depth (reads)

dim(DP) shows dimensions

quantile(DP)

Can type in console ?nameofpackage to learn about it

subset using [], [rows,columns] example: df[,c(1,4)]

one per line:
(x,
y,
z)

~/ goes to your home directory
```

### 09-16 Diversity Differentiation

We created a manhattan plot of the Fst values on 8 of the chromosomes within our Centaurea data. It shows where recombination happened the most in the dataset, and gives an overview of how to process a large set of chromosomes into a few specific ones.

%in% also found in ---\> looks for only the things in common (columns in this example)

```         
meta2 <- meta[meta$id %in% colnames(vcf@gt[,-1]),] 
# the [,] will let you choose what rows and columns to look at
-1 takes out the first column, which is just the format here
```

str() is another way to view the data, which gives you a lot more of the stats

```         
unique(vcf.div$CHROM) will give the unique values in the data
cbind binds columns from different sources
left_join joins 2 data frames
vcf.div.MHplot$V2 = as.numeric(vcf.div.MHplot$V2) makes the variable into its number version
suggestiveline = shows a cutoff line at a certain quantile of the data
```

### 09-24 Diversity Differentiation Cont.

Outliers on the Manhattan plot may be genes of interest potentially due to selection.

tidy operations use %\>% as a pipe for results into the next step

!= does not equal

& combines things (ex. filter(value!=0 & value \<0.25)

We created a long format of the Regions and then processed them into a histogram using ggplot to show Hs values. The graph shows how many SNPs are counted at a value of Hs, and the high number of 0s suggests no polymorphism at that site (only one allele present across all the samples)

### 09-26 PCA and Admixture

We continued our PCA analysis of the Centauria data, after we thinned the data.

A screeplot shows the magnitude of the eigenvalues in descending order, starting with the first PC value. We only want the first few PC values, which sum up most of the genetic variation

show() gives you an overview of your data

options(bitmapType = "cairo") helps when the plot doesn't work

### 10-1-2024 Admixture and Selection

We ran and plotted an admixture structure analysis of the Centaurea data to compare the fraction of ancestral DNA in individuals across different regions. We first compared the cross entropy values for different values of K to find one that best modeled our dataset, and then set that K value for our admixture to see the comparisons at that clustering.

After the selection lecture, we also modeled calculable selection outliers in R, looking at pvalues for PCAdapts on different chromosomes in a Manhattan plot.

entropy =T crossvalidates the K values

repetitions repeats the analysis a given number of times at a given value of K

par(mfrow= c(2,1)) stacks two plots so you can compare them

have to turn it off with dev.off()

Q scores here = fractional ancestral comparison for individuals

cbind pastes columns from two datasets next to eachother

min.maf is the minimum minor allele frequency
