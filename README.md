# RNA-Seq analysis pipeline

## Data
RNA-Seq data of patients who were diagnosed with acute leukemia and treated in Seoul St. Mary’s Hospital from February 2010 to March 2016. Standard diagnosis was established according to the WHO Classification of Tumours of Haematopoietic and Lymphoid Tissues based on bone marrow (BM) morphology, immunophenotyping, cytogenetic, and molecular genetic analysis.

> Ref: https://www.frontiersin.org/articles/10.3389/fonc.2021.717616/full


> RNAseq data: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA733693

Instrument: Illumina HiSeq 2500 - RNA-Seq

Reference genomes: Human GRCh38.p14

> https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/

## Set-up
- data/
    - SRR14689340_lymphoblastic.fastq
    - SRR14689341_lymphoblastic.fastq
    - SRR14689344_myeloid.fastq
    - SRR14689345_myeloid.fastq
    - GCF_000001405.40_GRCh38.p14_genomic.gtf
    - GCF_000001405.40_GRCh38.p14_genomic.fna
- run.sh
- deseq2.r
## Tools: 
sratoolskit, fastqc, multiqc, catadapt, bowtie2, samtools, tophat, deseq2

## Usage:
```
chmod +x run.sh

./run.sh
```

or submit run.sh to Job Scheduler
## Results:
1. Preprocessing: (run.sh) 
    * **1.multiqc_raw_data.pdf**: report for quality of RNA-Seq raw data
    * **2.multiqc_clean_data.pdf**: report for quality of RNA-Seq data after cleaning
    * **3.multiqc_tophat.pdf**: report for genes mapping
    * **4.multiqc_featuresCount.pdf**: report for genes counting
    * **features_matrix.txt**: gene counts for each sample
    * **conditions.txt**: condition of each sample
2. Differential gene expression analysis: (deseq2.r)
    * **5.MA_plot.pdf**: MA plot
    * **5.MA_plot_shrunken.pdf**: MA plot for shrunken log2 fold changes
    * **6.heatmap_count_ntd.pdf**: Heatmap of gene counts with Normalized Transformation
    * **6.heatmap_count_rlog.pdf**: Heatmap of gene counts with Regularized Logarithm Transformation
    * **6.heatmap_count_vst.pdf**: Heatmap of gene counts with Variance Stabilizing Transformation
    * **6.heatmap_sample_distance.pdf**: Heatmap of sample-to-sample distance
    * **7.PCA_samples.pdf**: Principle component plot of the samples
    * **8.gene_counts.pdf**: Plot counts of a specific gene
    * **lymph_myelo_results.csv**: Full report of differential expression analysis
    * **lymph_myelo_results_sig.csv**: Full report of differential expression analysis after apply pvalue cutoff of 0.05 and log2FoldChange of 4
    * **DE_genes.txt**: List of significant differential expression genes
3. Pathway analysis: (DAVID)
    * **9.DAVID_report.pdf**: Report from DAVID for the DE_genes.txt
    * **9.DAVID_KEGG_pathway.pdf**: Possible pathways from KEGG
    * **9.DAVID_reactome_pathway.pdf**: Possible pathways from reactome
