# scRNAalign_gene-dedup_count
Single Cell RNAseq pipeline readtransforming using umis, alignment, gene-deduplication by umi-tools and counting

pipeline was developed for CELSEQ-2 reads on the illumina platform but the use of UMIS package to transform reads should allow most single cell protocols to be run through this pipeline.

Read transformation using UMIS to combine reads into one containing cell, sample and umi sequences in the read name + an unique identifier (UID) created by concatenating sample cells and umi barcodes. This UID allows UMI-tools deduplication on bamfiles containing multiple cells.

The read name after transformation will include
CELL_BARCODE:UMI_BARCODE:SAMPLE_BARCODE:UID_[[samplebarcode][cellbarcode][umi]]

1 quality metrics from fastQC

2 UMIS read transformation
    -fastqtrasnform
    -cb_filter (filtering cellular barcoded on a predifined list)
    -sb_filter (filtering sample barcodes with edit distance 1 to have the same sample barcode in read name)
    -mb-filter (removing any umi reads that have a non ACGT base e.g. N)
    -add_uid (add the UID and save as fastq.gz)
    
3 Alignment using STAR

4 Read counting using featurecounts

5 adding XF:Z: tag to the BAM file containing the GENEID

6 deduplication using UMI-Tools

7 Generation of expression Matrix


There are options to change the analysis after aligment to count/dedup per contig instead of gene, or to skip deduplication

