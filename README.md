# RNA-Seq analysis pipeline

## Sample data
RNA-Seq data of patients who were diagnosed with acute leukemia and treated in Seoul St. Mary’s Hospital from February 2010 to March 2016. Standard diagnosis was established according to the WHO Classification of Tumours of Haematopoietic and Lymphoid Tissues based on bone marrow (BM) morphology, immunophenotyping, cytogenetic, and molecular genetic analysis.

Ref: https://www.frontiersin.org/articles/10.3389/fonc.2021.717616/full


RNAseq data: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA733693

Instrument: Illumina HiSeq 2500 - RNA-Seq

Reference genomes: Human GRCh38.p14

https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/

## Set-up
I have the following files in "/data" directory:

data/

-- SRR14689338_mixed.fastq

-- SRR14689339_mixed.fastq

-- SRR14689340_lymphoblastic.fastq

-- SRR14689341_lymphoblastic.fastq

-- SRR14689344_myeloid.fastq

-- SRR14689345_myeloid.fastq

-- GCF_000001405.40_GRCh38.p14_genomic.gtf

-- GCF_000001405.40_GRCh38.p14_genomic.fna


## Tools: 
fastqc, multiqc, bowtie2, tophat, deseq2

## Usage:
```
chmod +x run.sh

./run.sh
```

or submit run.sh to Job Scheduler