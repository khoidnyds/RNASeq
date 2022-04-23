#!/bin/sh
#SBATCH --time=10:00:00
#SBATCH --partition=normal_q
#SBATCH -n 32
#SBATCH --mem=100G
#SBATCH --account=aipmm

export PATH=/home/khoidnyds/RNAseq_old/tools/bowtie2-2.4.5:$PATH
export PATH=/home/khoidnyds/RNAseq_old/tools/tophat-2.1.1:$PATH

export THREADS=32
# fasterq-dump --outdir data --mem 10G --split-3 --threads $THREADS --skip-technical  --print-read-nr SRR14689338 SRR14689339 SRR14689340 SRR14689341 SRR14689344 SRR14689345

export REF="data/GCF_000001405.40_GRCh38.p14_genomic.fna"
export REF_BUILD="5.bowtie2/REF_GENOME"
export ANNO="data/GCF_000001405.40_GRCh38.p14_genomic.gtf"

export MIXED_1="data/SRR14689338_mixed.fastq"
export MIXED_2="data/SRR14689339_mixed.fastq"
export LYMPHOBLASTIC_1="data/SRR14689340_lymphoblastic.fastq"
export LYMPHOBLASTIC_2="data/SRR14689341_lymphoblastic.fastq"
export MYELOID_1="data/SRR14689344_myeloid.fastq"
export MYELOID_2="data/SRR14689345_myeloid.fastq"

export MIXED_1_CLEAN="3.cutadapt/SRR14689338_mixed_clean.fastq"
export MIXED_2_CLEAN="3.cutadapt/SRR14689339_mixed_clean.fastq"
export LYMPHOBLASTIC_1_CLEAN="3.cutadapt/SRR14689340_lymphoblastic_clean.fastq"
export LYMPHOBLASTIC_2_CLEAN="3.cutadapt/SRR14689341_lymphoblastic_clean.fastq"
export MYELOID_1_CLEAN="3.cutadapt/SRR14689344_myeloid_clean.fastq"
export MYELOID_2_CLEAN="3.cutadapt/SRR14689345_myeloid_clean.fastq"

export MIXED_1_SORT="7.samtools/SRR14689338_mixed.bam"
export MIXED_2_SORT="7.samtools/SRR14689339_mixed.bam"
export LYMPHOBLASTIC_1_SORT="7.samtools/SRR14689340_lymphoblastic.bam"
export LYMPHOBLASTIC_2_SORT="7.samtools/SRR14689341_lymphoblastic.bam"
export MYELOID_1_SORT="7.samtools/SRR14689344_myeloid.bam"
export MYELOID_2_SORT="7.samtools/SRR14689345_myeloid.bam"

export MIXED_1_FC="8.features/mixed_1.txt"
export MIXED_2_FC="8.features/mixed_2.txt"
export LYMPHOBLASTIC_1_FC="8.features/lymphoblastic_1.txt"
export LYMPHOBLASTIC_2_FC="8.features/lymphoblastic_2.txt"
export MYELOID_1_FC="8.features/myeloid_1.txt"
export MYELOID_2_FC="8.features/myeloid_2.txt"

# STEP 1. QUALITY CONTROL: remove reads being shorter than 20 nucleotides, reads having quality score smaller than 20
# mkdir 1.fastqc
# fastqc -t 32 data/*.fastq* -o 1.fastqc
# multiqc 1.fastqc -f -o 2.multiqc

# mkdir 3.cutadapt
# cutadapt --cores=0 --minimum-length 20 --quality-cutoff 20 -o $MIXED_1_CLEAN $MIXED_1 > 3.cutadapt/report_mix_1.txt
# cutadapt --cores=0 --minimum-length 20 --quality-cutoff 20 -o $MIXED_2_CLEAN $MIXED_2 > 3.cutadapt/report_mix_2.txt
# cutadapt --cores=0 --minimum-length 20 --quality-cutoff 20 -o $LYMPHOBLASTIC_1_CLEAN $LYMPHOBLASTIC_1 > 3.cutadapt/report_lymphoblastic_1.txt
# cutadapt --cores=0 --minimum-length 20 --quality-cutoff 20 -o $LYMPHOBLASTIC_2_CLEAN $LYMPHOBLASTIC_2 > 3.cutadapt/report_lymphoblastic_2.txt
# cutadapt --cores=0 --minimum-length 20 --quality-cutoff 20 -o $MYELOID_1_CLEAN $MYELOID_1 > 3.cutadapt/report_myeloid_1.txt
# cutadapt --cores=0 --minimum-length 20 --quality-cutoff 20 -o $MYELOID_2_CLEAN $MYELOID_2 > 3.cutadapt/report_myeloid_2.txt
# multiqc 3.cutadapt -f -o 4.multiqc


# STEP 2. MAPPING: minimum length of segments of read is 18, disable coverage search for junctions,
# mkdir 5.bowtie2
# bowtie2-build $REF $REF_BUILD --threads $THREADS
# cp $REF 5.bowtie2
# mv 5.bowtie2/GCF_000001405.40_GRCh38.p14_genomic.fna 5.bowtie2/REF_GENOME.fa

mkdir 6.tophat2
tophat2 -o 6.tophat2/mix_1 --num-threads $THREADS --segment-length 18 --no-coverage-search --mate-inner-dist 200 $REF_BUILD $MIXED_1_CLEAN
tophat2 -o 6.tophat2/mix_2 --num-threads $THREADS --segment-length 18 --no-coverage-search --mate-inner-dist 200 $REF_BUILD $MIXED_2_CLEAN
tophat2 -o 6.tophat2/lymphoblastic_1 --num-threads $THREADS --segment-length 18 --no-coverage-search --mate-inner-dist 200 $REF_BUILD $LYMPHOBLASTIC_1_CLEAN
tophat2 -o 6.tophat2/lymphoblastic_2 --num-threads $THREADS --segment-length 18 --no-coverage-search --mate-inner-dist 200 $REF_BUILD $LYMPHOBLASTIC_2_CLEAN
tophat2 -o 6.tophat2/myeloid_1 --num-threads $THREADS --segment-length 18 --no-coverage-search --mate-inner-dist 200 $REF_BUILD $MYELOID_1_CLEAN
tophat2 -o 6.tophat2/myeloid_2 --num-threads $THREADS --segment-length 18 --no-coverage-search --mate-inner-dist 200 $REF_BUILD $MYELOID_2_CLEAN


# STEP 3. SORTING
# mkdir 7.samtools
# samtools sort -n -o $MIXED_1_SORT -@ $THREADS 6.tophat2/mix_1/accepted_hits.bam
# samtools index $MIXED_1_SORT
# samtools sort -n -o $MIXED_2_SORT -@ $THREADS 6.tophat2/mix_2/accepted_hits.bam
# samtools index $MIXED_2_SORT
# samtools sort -n -o $LYMPHOBLASTIC_1_SORT -@ $THREADS 6.tophat2/lymphoblastic_1/accepted_hits.bam
# samtools index $LYMPHOBLASTIC_1_SORT
# samtools sort -n -o $LYMPHOBLASTIC_2_SORT -@ $THREADS 6.tophat2/lymphoblastic_2/accepted_hits.bam
# samtools index $LYMPHOBLASTIC_2_SORT
# samtools sort -n -o $MYELOID_1_SORT -@ $THREADS 6.tophat2/myeloid_1/accepted_hits.bam
# samtools index $MYELOID_1_SORT
# samtools sort -n -o $MYELOID_2_SORT -@ $THREADS 6.tophat2/myeloid_2/accepted_hits.bam
# samtools index $MYELOID_2_SORT

# STEP 4. COUNTING: Count reads for each gene based on the sorted bam files
# mkdir 8.features
# featureCounts -t exon -g gene_id -T $THREADS -a $ANNO -o $MIXED_1_FC $MIXED_1_SORT
# featureCounts -t exon -g gene_id -T $THREADS -a $ANNO -o $MIXED_2_FC $MIXED_2_SORT
# featureCounts -t exon -g gene_id -T $THREADS -a $ANNO -o $LYMPHOBLASTIC_1_FC $LYMPHOBLASTIC_1_SORT
# featureCounts -t exon -g gene_id -T $THREADS -a $ANNO -o $LYMPHOBLASTIC_2_FC $LYMPHOBLASTIC_2_SORT
# featureCounts -t exon -g gene_id -T $THREADS -a $ANNO -o $MYELOID_1_FC $MYELOID_1_SORT
# featureCounts -t exon -g gene_id -T $THREADS -a $ANNO -o $MYELOID_2_FC $MYELOID_2_SORT


# STEP 5. Differential expression analysis (DE) can be done by using a R package: DESeq2. (https://www.r-project.org/, https://www.rstudio.com/products/rstudio/, and https://www.bioconductor.org/packages/release/bioc/html/DESeq2.html). Try to describe the codes from the DESeq2 to generate different expressed gene list, MA-plot, PCA plot, heatmap for the DE genes and ggplot to display any gene expression level (DESeq2 vignette will be helpful).
# Differential expression genes
# mkdir 6.countReadPerGene
# htseq-count -t exon -i gene_id -f bam $F77_BAM $F78_BAM $F80_BAM $F81_BAM $ANNO > 6.countReadPerGene/readCount_raw.txt
# multiqc 6.countReadPerGene -f -o 7.multiqc
# head -n -5 6.countReadPerGene/readCount_raw.txt > 6.countReadPerGene/readCount.txt
# Rscript deseq2.r

# STEP 6. Use web-based pathway analysis (https://david.ncifcrf.gov/tools.jsp) and find possible biological pathways that related to the DE genes.