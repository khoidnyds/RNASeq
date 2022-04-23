#!/bin/sh
#SBATCH --time=10:00:00
#SBATCH --partition=normal_q
#SBATCH -n 32
#SBATCH --mem=100G
#SBATCH --account=aipmm

export PATH=/home/khoidnyds/RNAseq_old/tools/bowtie2-2.4.5:$PATH
export PATH=/home/khoidnyds/RNAseq_old/tools/tophat-2.1.1:$PATH

export THREADS=32

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

export MIXED_1_SORT="8.samtools/SRR14689338_mixed.bam"
export MIXED_2_SORT="8.samtools/SRR14689339_mixed.bam"
export LYMPHOBLASTIC_1_SORT="8.samtools/SRR14689340_lymphoblastic.bam"
export LYMPHOBLASTIC_2_SORT="8.samtools/SRR14689341_lymphoblastic.bam"
export MYELOID_1_SORT="8.samtools/SRR14689344_myeloid.bam"
export MYELOID_2_SORT="8.samtools/SRR14689345_myeloid.bam"

export FEATURES_COUNT="9.features/features_count.txt"
export FEATURES_MATRIX="9.features/features_matrix.txt"

# STEP 0: GET DATA: from SRA
fasterq-dump --outdir data --mem 10G --split-3 --threads $THREADS --skip-technical --print-read-nr SRR14689338 SRR14689339 SRR14689340 SRR14689341 SRR14689344 SRR14689345


# STEP 1. QUALITY CONTROL: remove reads being shorter than 20 nucleotides, reads having quality score smaller than 20
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


# STEP 2. MAPPING: minimum length of segments of read is 18, disable coverage search for junctions,
mkdir 5.bowtie2
bowtie2-build $REF $REF_BUILD --threads $THREADS
cp $REF 5.bowtie2
mv 5.bowtie2/GCF_000001405.40_GRCh38.p14_genomic.fna 5.bowtie2/REF_GENOME.fa

mkdir 6.tophat2
tophat2 -o 6.tophat2/mix_1 --num-threads $THREADS --segment-length 18 --no-coverage-search --mate-inner-dist 200 $REF_BUILD $MIXED_1_CLEAN
tophat2 -o 6.tophat2/mix_2 --num-threads $THREADS --segment-length 18 --no-coverage-search --mate-inner-dist 200 $REF_BUILD $MIXED_2_CLEAN
tophat2 -o 6.tophat2/lymphoblastic_1 --num-threads $THREADS --segment-length 18 --no-coverage-search --mate-inner-dist 200 $REF_BUILD $LYMPHOBLASTIC_1_CLEAN
tophat2 -o 6.tophat2/lymphoblastic_2 --num-threads $THREADS --segment-length 18 --no-coverage-search --mate-inner-dist 200 $REF_BUILD $LYMPHOBLASTIC_2_CLEAN
tophat2 -o 6.tophat2/myeloid_1 --num-threads $THREADS --segment-length 18 --no-coverage-search --mate-inner-dist 200 $REF_BUILD $MYELOID_1_CLEAN
tophat2 -o 6.tophat2/myeloid_2 --num-threads $THREADS --segment-length 18 --no-coverage-search --mate-inner-dist 200 $REF_BUILD $MYELOID_2_CLEAN
multiqc 6.tophat2 -f -o 7.multiqc


# STEP 3. SORTING: Sort alignment file by read name.
mkdir 8.samtools
samtools sort -@ $THREADS -n -o $MIXED_1_SORT  6.tophat2/mix_1/accepted_hits.bam
samtools sort -@ $THREADS -n -o $MIXED_2_SORT  6.tophat2/mix_2/accepted_hits.bam
samtools sort -@ $THREADS -n -o $LYMPHOBLASTIC_1_SORT 6.tophat2/lymphoblastic_1/accepted_hits.bam
samtools sort -@ $THREADS -n -o $LYMPHOBLASTIC_2_SORT 6.tophat2/lymphoblastic_2/accepted_hits.bam
samtools sort -@ $THREADS -n -o $MYELOID_1_SORT 6.tophat2/myeloid_1/accepted_hits.bam
samtools sort -@ $THREADS -n -o $MYELOID_2_SORT 6.tophat2/myeloid_2/accepted_hits.bam


# STEP 4. COUNTING: Count reads for each gene based on the sorted bam files, minimum mapping quality is 10,
mkdir 9.features
featureCounts -Q 10 -t exon -g gene_id -T $THREADS -a $ANNO -o $FEATURES_COUNT $MIXED_1_SORT $MIXED_2_SORT $LYMPHOBLASTIC_1_SORT $LYMPHOBLASTIC_2_SORT $MYELOID_1_SORT $MYELOID_2_SORT
awk -F'\t' '{ print $1,"\t",$7,"\t",$8,"\t",$9,"\t",$10,"\t",$11,"\t",$12 }' $FEATURES_COUNT > $FEATURES_MATRIX
multiqc 9.features -f -o 10.multiqc


# STEP 5. DIFFERENTIAL EXPRESSION ANALYSIS: the input is features_matrix.txt and conditions.txt
Rscript deseq2.r


# STEP 6. PATHWAY ANALYSIS (https://david.ncifcrf.gov/tools.jsp) and find possible biological pathways that related to the DE genes.
# Attached images