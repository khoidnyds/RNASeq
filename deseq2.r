# Try to describe the codes from the DESeq2 to generate 
# different expressed gene list, 
# MA-plot
# PCA plot, 
# heatmap for the DE genes and 
# ggplot to display any gene expression level (DESeq2 vignette will be helpful)


library(DESeq2)
library("pheatmap")
library("RColorBrewer")
library(ggplot2)

pvalue_cutoff = 0.05
log2FoldChange_cutoff = 4

# read features_matrix (49677 x 4) 
cts <- read.csv("features_matrix.txt", sep = '\t', header=T,  row.names = "Geneid")
cts <- as.matrix(cts)

# read conditions (4 x 1)
coldata <- read.csv("conditions.txt", row.names = 1)
coldata$condition <- factor(coldata$condition)

# create DESeqDataSet obj 
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)
# run DE 
dds <- DESeq(dds)

res <- results(dds) # 49677 rows and 6
resultsNames(dds)
resLFC <- lfcShrink(dds, coef="condition_myelo_vs_lymph", type='apeglm')
# remove row if pvalue is N/A
res <- na.omit(res) # 17209 rows and 6 columns 

# plot MA, highlight the dot having pvalue smaller than the cutoff
par(mar = c(5, 5, 5, 5))
plotMA(res, alpha=pvalue_cutoff, main="MA plot. Many genes (blue dots) have the pvalue < 0.05. \nThe number of upregulated genes are similar to the number of downregulated genes.")
plotMA(resLFC, alpha=pvalue_cutoff, , main="MA plot for the shrunken log2 fold changes. Many genes (blue dots) have the pvalue < 0.05. \nThe number of upregulated genes are similar to the number of downregulated genes.")

# order the results based on pvalue in ascending order
resOrdered <- res[order(res$pvalue),]
write.csv(as.data.frame(resOrdered), file="lymph_myelo_results.csv") # 17209 genes

# select only significant gene based on pvalue threshold and log2FoldChange threshold
resSig <- res[res$padj < pvalue_cutoff, ]
resSig['abs_log2FC'] <- abs(resSig$log2FoldChange)
resSig <- resSig[resSig$abs_log2FC > log2FoldChange_cutoff, ]
resSig <- resSig[order(resSig$pvalue),]
write.csv(as.data.frame(resSig), file="lymph_myelo_results_sig.csv")
write(rownames(resSig), "DE_genes.txt", sep = "\n") # 595 genes

# Visualize heatmap of the count matrix of most 20 common genes
ntd <- normTransform(dds)
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:30]
df <- as.data.frame(colData(dds)[,c("condition")])
colnames(df) <- c("Lympho vs Myelo")
rownames(df) <- c("Lymph_1","Lymph_2","Myelo_1","Myelo_2")
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=FALSE, annotation_col=df,
         main="Heatmap of count matrix (Normalized Transformation) of most 20 common genes.\nSome genes show the signification of gene expression between lympho and myelo groups\nsuch as AHNAK, XIST, CHST15, VCAN, LYST.")

vsd <- vst(dds, blind=FALSE)
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=FALSE, annotation_col=df,
         main="Heatmap of count matrix (Variance Stabilizing Transformation) of most 20 common genes.\nSome genes show the signification of gene expression between lympho and myelo groups\nsuch as AHNAK, SORL1, XIST, CHST15, VCAN, LYST.")

rld <- rlog(dds, blind=FALSE)
pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=FALSE, annotation_col=df,
         main="Heatmap of count matrix (Regularized Logarithm Transformation) of most 20 common genes.\nSome genes show the signification of gene expression between lympho and myelo groups\nsuch as AHNAK, XIST, CHST15, VCAN, LYST.")

# Visualize heatmap of the sample-to-sample distances
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors,
         main="Heatmap of sample-to-sample distances. Distances between samples\nin a same group are smaller than distances between samples in different groups.")

# PCA of samples
pcaData <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=2) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("PCA of the samples. PC1 explains 60% variance of the data.\nPC2 explains 29% variance of the data. The Lympho samples are in one cluster.")
  coord_fixed() 

# Display gene expression level of any gene
plotCounts(dds, gene="POU4F1 ", intgroup="condition")
