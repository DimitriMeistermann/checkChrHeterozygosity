# Description
This pipeline aims to provide a simple alignment tool for  RNAseq data (without UMI). The fastq files are aligned on a reference genome with STAR and counted with HTSeq-Count. 

# Prerequisites
- Git
- Conda


# Input
The *fastq.gz* files need to be gathered in a directory. The pathway to this directory will be specified in the config.json. **One fastq** file per sample (or two if data are sequenced in paired end).


# Installation
~~~
git git@github.com:DimitriMeistermann/simpleSNPcall.git
cd simpleSNPcall
conda env create -f virtualEnvs/simpleCallSNP.yml
~~~

# Configuration

- **IS_PAIRED_END**: *true* or *false*, is the data in paired end (forward and reverse reads ) ?
- **THREAD_PER_SAMPLE**: Number of cores used by most of the workflow jobs.
- **FASTA_REFERENCE**: Path of fasta file of the  reference genome
- **GTF_REFERENCE**: Path of the GTF file (feature of the reference genome)
- **FASTQ_PATH**: Directory that contains the input fastqs.
- **OUTPUT_PATH**: Path where the results will be written.
- **PAIR_END_FILE_PATTERN**: characters between sample name and 1 or 2 for paired end data. Example: if forward reads are in *sample_1.fastq.gz* and reverse in *sample_2.fastq.gz*, then *PAIR_END_FILE_PATTERN = "_"*
- **USE_TRIMMOMATIC**  *true* or *false*, use trimmomatic for removing adpataters/bad quality part from the FASTQs
- **FEATURE_TYPE**: Feature of the GTF (second column) were aligned reads are counted (example: Exon, Gene,...).
- **FEATURE_ID**: ID that will be the rows of the count table (gene_name for gene symbols).

# Execution of the pipeline
## Normal launch (with 8 cores), with default config.json
~~~
conda activate simpleSNP
snakemake -rp -j 8
~~~

## With specific config file
~~~
conda activate simpleSNP
snakemake -rp -j 8 --configfile configFileSave/[yourConfig.json]
~~~

# Output data
- **log**: log files from executed jobs, in the form of *RULE_NAME[_sampleName].(err|out)*. Additional logs are maybe saved depending the tool that was used.
- **FASTQ_TRIM**: optional, if trimmomatic is used, trimmed fastqs are stored in this folder.
- **FASTQ_TRIM_UNPAIR**: if trimmomatic is used with paired end FASTQs, reads that lost their pair are stored here.
- **fastQC**: Zip and HTML results from the execution of FASTQC.
- **BAM**: this directory contains the results of alignment in the form of BAM files.
- **SORTED_BAM**: sorted/indexed BAM'
- **VCF**: Result of variant calling in VCF format
- **TSV**: Result of variant calling in TSV format
- **counts**: this directory contains the resulting counts files obtained from htseqcount.
- **results**: this directory contains the raw counts table, a table of alignment stats, the multiqc reports that present different quality scores from the data, the merged TSV file of variant calling results, and a barplot of allelic distribution of variants.
