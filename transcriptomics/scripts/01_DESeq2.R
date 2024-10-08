### Code for analyzing RNAseq (Transcriptomics) data using DEseq2

# load libraries
library(DESeq2)
library(ggplot2)

setwd("~/projects/eco_genomics/transcriptomics/")

# Import counts matrix

countsTable <- read.table("/gpfs1/cl/pbio3990/Transcriptomics/tonsa_counts.txt", 
                          header = TRUE, row.names = 1)

countsTableRound <- round(countsTable) # because DEseq2 doesn't like decimals
tail(countsTableRound) #see bottom of data

conds <- read.delim("/gpfs1/cl/pbio3990/Transcriptomics/experimental_details.txt",
                    header = TRUE, stringsAsFactors = TRUE, row.names = 1)
head(conds)

dds <- DESeqDataSetFromMatrix(countData= countsTableRound, colData = conds, 
                              design = ~ DevTemp + FinalTemp)
