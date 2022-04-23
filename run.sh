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

export REFERENCE="data/GCF_000001405.40_GRCh38.p14_genomic.fna"
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

export ALIGNMENT_SORTED_77="3.samtools/alignment_sorted_77.bam"
export ALIGNMENT_SORTED_80="3.samtools/alignment_sorted_80.bam"
export FEATURES_COUNT_77="4.features/count_77.txt"
export FEATURES_COUNT_80="4.features/count_80.txt"

# STEP 1. QUALITY CONTROL: remove reads being shorter than 20 nucleotides, read having quality score smaller than 20
mkdir 1.fastqc
fastqc -t 32 data/*.fastq* -o 1.fastqc
multiqc 1.fastqc -f -o 2.multiqc
mkdir 3.cutadapt
cutadapt --cores=0 --minimum-length 20 --quality-cutoff 20 -o $MIXED_1_CLEAN $MIXED_1 > 3.cutadapt/report_mix_1.txt
cutadapt --cores=0 --minimum-length 20 --quality-cutoff 20 -o $MIXED_2_CLEAN $MIXED_2 > 3.cutadapt/report_mix_2.txt
cutadapt --cores=0 --minimum-length 20 --quality-cutoff 20 -o $LYMPHOBLASTIC_1_CLEAN $LYMPHOBLASTIC_1 > 3.cutadapt/report_lymphoblastic_1.txt
cutadapt --cores=0 --minimum-length 20 --quality-cutoff 20 -o $LYMPHOBLASTIC_2_CLEAN $LYMPHOBLASTIC_2 > 3.cutadapt/report_lymphoblastic_2.txt
cutadapt --cores=0 --minimum-length 20 --quality-cutoff 20 -o $MYELOID_1_CLEAN $MYELOID_1 > 3.cutadapt/report_myeloid_1.txt
cutadapt --cores=0 --minimum-length 20 --quality-cutoff 20 -o $MYELOID_2_CLEAN $MYELOID_2 > 3.cutadapt/report_myeloid_2.txt
multiqc 3.cutadapt -f -o 4.multiqc


# STEP 2. MAPPING
mkdir 5.bowtie2
bowtie2-build $REFERENCE 5.bowtie2/REF_GENOME --threads $THREADS
cp $REFERENCE 5.bowtie2/REF_GENOME.fa
mkdir 6.tophat2
# tophat2 -o 2.tophat2_77 --num-threads $THREADS --segment-length 18 --no-coverage-search --mate-inner-dist 200 2.bowtie2/DM_genome $F77_1 $F77_2
# tophat2 -o 2.tophat2_80 --num-threads $THREADS --segment-length 18 --no-coverage-search --mate-inner-dist 200 2.bowtie2/DM_genome $F80_1 $F80_2


# # # STEP 3. Use samtools to sort the alignment files (bam format) by read name. (http://www.htslib.org/doc/samtools.html)
# # mkdir 3.samtools
# samtools sort -n -o $ALIGNMENT_SORTED_77 -@ $THREADS 2.tophat2_77/accepted_hits.bam
# samtools sort -n -o $ALIGNMENT_SORTED_80 -@ $THREADS 2.tophat2_80/accepted_hits.bam


# # # STEP 4. Count reads for each gene based on the sorted bam files. (https://htseq.readthedocs.io/en/release_0.11.1/count.html)
# # mkdir 4.features
# featureCounts -Q 10 -t exon -g gene_id -a $ANNO -o $FEATURES_COUNT_77 $ALIGNMENT_SORTED_77
# featureCounts -Q 10 -t exon -g gene_id -a $ANNO -o $FEATURES_COUNT_80 $ALIGNMENT_SORTED_80


# # # STEP 5. Differential expression analysis (DE) can be done by using a R package: DESeq2. (https://www.r-project.org/, https://www.rstudio.com/products/rstudio/, and https://www.bioconductor.org/packages/release/bioc/html/DESeq2.html). Try to describe the codes from the DESeq2 to generate different expressed gene list, MA-plot, PCA plot, heatmap for the DE genes and ggplot to display any gene expression level (DESeq2 vignette will be helpful).


# # # STEP 6. Use web-based pathway analysis (https://david.ncifcrf.gov/tools.jsp) and find possible biological pathways that related to the DE genes.