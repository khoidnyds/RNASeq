# Try to describe the codes from the DESeq2 to generate different expressed gene list, 
# MA-plot
# PCA plot, 
# heatmap for the DE genes and 
# ggplot to display any gene expression level (DESeq2 vignette will be helpful)


# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("airway")
library(DESeq2)
# read from 12.coun 
cts <- read.csv("6.countReadPerGene/readCount.txt", sep = '\t', header = F, row.names = 1)
cts <- as.matrix(cts)

coldata <- read.csv("6.countReadPerGene/conditions.txt", row.names = 1)
coldata$condition <- factor(coldata$condition)

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)
# Differential expression analysis
dds <- DESeq(dds)
res <- results(dds, contrast=c("condition","treated","untreated"))

#pdf("MA_plot.pdf")
fig1 <- plotMA(res, ylim=c(-2,2))
png("fig1.png", width = 500, height = 500, res = 300)
plot(fig1, xlim=c(0,100), ylim=c(-2,2))
dev.off() 


# res_filter <- res[res$pvalue < 0.05, ]
# res_filter <- na.omit(res_filter)
# res_filter <- res_filter[res_filter$abs_log2FC > 1, ]


# write(rownames(res_filter), "8.DE/DE.txt", sep = "\n")


# library("pheatmap")
# select <- order(rowMeans(counts(res,normalized=TRUE)),
#                 decreasing=TRUE)[1:20]
# df <- as.data.frame(colData(dds)[,c("condition","type")])
# pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
#          cluster_cols=FALSE, annotation_col=df)