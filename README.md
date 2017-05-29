# _scRNApipe_

The scRNApipe pipeline was originally designed to preprocess and analyse scRNA-Seq data, following the CEL-Seq2 protocol, on the Illumina platforms. Nevertheless, by data transformation, provided by UMIS package, will allow most of single cell protocols to be run through this pipeline.


Read transformation will combine reads into one containing cell, sample and umi barcode sequences incorporated in the read name + a unique identifier (UID) created by concatenating those three barcodes. This UID will allow UMI-tools to remove PCR duplicates on bamfiles containing multiple cells.

In principle, the raw data will be readtransformed and filtered (by UMIs), aligned, gene-deduplicate (by UMI-tools) and counted.

 * Quality metrics (optional)
 * Preprocessing 
 * Aligning
 * Main Analysis
 * Expression Matrix

>Main _options_ to tune in this pipeline:
>	* The Main Analysis can run in count/dedup per contig or default mode (instead of gene)
>* Skip deduplication

<br />

## scRNApipe In Details

#### 1. Quality metrics from FastQC
Detailed reports will be generated for each sample by FastQC. An the end a summarised report will be available for an overall review of all samples at once. 

#### Preprocessing the reads


2.	_umis fastqtrasnform_ (read transformation)
3.	_cb_filter_ (filtering reads with non-matching CELLULAR barcodes (<span style="color:red">CB</span>) | 1 mismatch is allowed) 
4.	_sb_filter_ (filtering reads with non-matching SAMPLE barcodes (<span style="color:blue">SB</span>) | 1 mismatch is allowed) 
5.  _mb_filter_ (removing reads with ambiguous (e.g N) bases in the UMI barcodes) 
6.  _add_uid_ (add the UID and save as fastq.gz)

The read name after preprocessing will include
CELL_BARCODE:UMI_BARCODE:SAMPLE_BARCODE:UID_[[<span style="color:blue">samplebarcode</span>][<span style="color:red">cellbarcode</span>][umi]]

#### 7. Alignment using STAR aligner
Aligning the preprocessed reads against the reference genome by the use of the STAR aligner.

#### Main Analysis

8.	Counting reads using featureCounts
9.	Adding XF:Z: tag to the BAM file containing the GeneID
10.	Deduplication using UMI-Tools

#### Generation of Expression Matrix
Generate the Expression Matrix based on the GeneID tags 

<br />

## __Installation__ and __Info__

If you'd like to work directly from the git repository:

	$ git clone https://github.com/MarinusVL ...

Enter repository and run:
	
    $ python setup.py install

### __Executing__

After installation the pipeline can be used:

	$ scRNApipe <configuration_file.txt>

### __Help__

For further information about each compartment of the pipeline you can run:

	$ scRNApipe --help

### __Dependencies__

scRNApipe is dependent on umis, umi_tools, numpy, pysam, STAR, featureCounts, fastqc and multiqc
and Python 2.7
