import glob, os, sys, json, shutil


WORKING_DIR = os.path.dirname(workflow.snakefile)
if config=={}:
	print("Default config file loaded, from " + WORKING_DIR + "/config.json")
	configfile: WORKING_DIR+"/config.json"

## creation of the logs subdirectory
if not os.path.exists(WORKING_DIR+"/log"):
	os.mkdir(WORKING_DIR+"/log")

#put all config variable as variable in the snakefile
for configVar in config:
	if isinstance(config[configVar], str): exec(configVar+"= '"+config[configVar]+"'")
	else: exec(configVar+"="+str(config[configVar]))


## test of the path provided in the config.json file
if not os.path.exists(FASTQ_PATH):
	print("The directory " + FASTQ_PATH + " doesn't exist. Check the field FASTQ_PATH into the config.json file.")
	sys.exit(0)
else:
	## If the path ends by /, the / is suppressed
	if ( FASTQ_PATH[-1:] == "/" ):
		FASTQ_PATH = FASTQ_PATH[:-1]

INPUT_FASTQS = glob.glob(FASTQ_PATH+'/*.fastq.gz')

SAMPLESsplitted = [os.path.basename(f).split(".") for f in INPUT_FASTQS]
SAMPLES=[]

#remove .fastq.gz to get sample names
for s in SAMPLESsplitted:
	SAMPLES.append(".".join(s[0:-2]))

if(OUTPUT_PATH[-1] == "/") : OUTPUT_PATH = OUTPUT_PATH[:-1]

## suppress the .R1. and .R2. elements for paired-end fastq files for the alignement processus in SAMPLES
if IS_PAIRED_END:
	SAMPLES = [itemR2 for itemR2 in SAMPLES if (PAIR_END_FILE_PATTERN+"2") not in itemR2]	
	SAMPLES = [itemR1.replace((PAIR_END_FILE_PATTERN+"1"),'') for itemR1 in SAMPLES]


if USE_TRIMMOMATIC:
	ALIGN_FASTQ_FOLDER = OUTPUT_PATH + "/FASTQ_TRIM"
else:
	ALIGN_FASTQ_FOLDER = FASTQ_PATH

if IS_PAIRED_END: PAIR_SUFFIX = [PAIR_END_FILE_PATTERN+"1",PAIR_END_FILE_PATTERN+"2"]
else: PAIR_SUFFIX = [""]


if IS_PAIRED_END:  print("Workflow set on paired end mode")
else : print("Workflow set on single end mode")


##############
rule all: 
	input: OUTPUT_PATH+"/results/multiqc_report.html"

if IS_PAIRED_END : 
	rule TRIMMOMATIC:
		input: expand(FASTQ_PATH+"/{{sample}}{pair}.fastq.gz",pair=PAIR_SUFFIX)
		output:
			pairedR1=OUTPUT_PATH+"/FASTQ_TRIM/{sample}"+PAIR_SUFFIX[0]+".fastq.gz",
			pairedR2=OUTPUT_PATH+"/FASTQ_TRIM/{sample}"+PAIR_SUFFIX[1]+".fastq.gz",
			unpairedR1=OUTPUT_PATH+"/FASTQ_TRIM_UNPAIR/{sample}"+PAIR_SUFFIX[0]+".fastq.gz",
			unpairedR2=OUTPUT_PATH+"/FASTQ_TRIM_UNPAIR/{sample}"+PAIR_SUFFIX[1]+".fastq.gz",
		shell: """
		trimmomatic PE {input} {output.pairedR1} {output.unpairedR1} {output.pairedR2} {output.unpairedR2} \
		ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36
		"""	

rule FASTQC:
	input: ALIGN_FASTQ_FOLDER+"/{sample}{pair}.fastq.gz"
	output: multiext(OUTPUT_PATH+"/fastQC/{sample}{pair}_fastqc",".zip",".html")
	shell: """
	fastqc -o {OUTPUT_PATH}/fastQC {input}
	"""
rule BWA_INDEX:
	input: FASTA_REFERENCE
	output: FASTA_REFERENCE+".bwt"
	shell: """
	bwa index {input}
	"""

# rule CreateSequenceDictionary:
# 	input: FASTA_REFERENCE
# 	output: FASTA_REFERENCE+".dict"
# 	shell: """
# 	picard CreateSequenceDictionary R={input} O={output}
# 	"""
### ALIGN EACH PAIR OF FASTQ FILES
rule ALIGN:
	input:
		fastq = expand(ALIGN_FASTQ_FOLDER+"/{{sample}}{pair}.fastq.gz",pair=PAIR_SUFFIX),
		index = FASTA_REFERENCE+".bwt"
	log: out=OUTPUT_PATH+"/log/ALIGN_{sample}.out"
	output:	OUTPUT_PATH+"/BAM/{sample}.bam"
	shell:	"""
	bwa mem -t 1 -M \
	-R '@RG\\tID:{wildcards.sample}\\tLB:{wildcards.sample}\\tSM:{wildcards.sample}\\tPL:Illumina\\tCN:CENTER' \
	{FASTA_REFERENCE} {input.fastq} 2> {log.out} | samtools view -bS > {output}
	"""		

rule SORT_BAM:
	input: OUTPUT_PATH+"/BAM/{sample}.bam"
	output: OUTPUT_PATH+"/SORTED_BAM/{sample}.bam"
	params: cpu = THREAD_PER_SAMPLE
	shell: """
	samtools sort -@{params.cpu} -o {output} {input} 
	"""

rule INDEX_BAM:
	input: OUTPUT_PATH+"/SORTED_BAM/{sample}.bam"
	output: OUTPUT_PATH+"/SORTED_BAM/{sample}.bam.bai"
	params: cpu = THREAD_PER_SAMPLE
	shell: """
	samtools index -@{params.cpu} {input}
	"""
	
###############
rule mpileUp:
	input:
		sorted = OUTPUT_PATH+"/SORTED_BAM/{sample}.bam",
		index = OUTPUT_PATH+"/SORTED_BAM/{sample}.bam.bai"
	output: OUTPUT_PATH+"/VCF/{sample}.vcf"
	shell: """
	bcftools mpileup -Ou -f {FASTA_REFERENCE} {input.sorted} | bcftools call -Ou -mv  | bcftools norm -Ou -f {FASTA_REFERENCE}  -o {output}
	"""

### FILTER SNVS
#rule filterSNVs:
#	input:
#		vcf = OUTPUT_PATH+"/VCF/{sample}.vcf",
#		dict = FASTA_REFERENCE+".dict"
#	output:	vcf = OUTPUT_PATH+"/filteredVCF/{sample}.vcf"
#	shell: """
#	gatk -Xmx3g -T VariantFiltration \
#	-R {FASTA_REFERENCE} \
#	--logging_level FATAL  \
#	--filterExpression 'QD < 2.0' --filterName 'FAILS_HARD_FILTER_SNP_QD' \
#	--filterExpression 'FS > 60.0' --filterName 'FAILS_HARD_FILTER_SNP_FS' \
#	--filterExpression 'MQ < 30.0' --filterName 'FAILS_HARD_FILTER_SNP_MQ' \
#	--filterExpression 'MQRankSum < -12.5' --filterName 'FAILS_HARD_FILTER_SNP_MQRankSum' \
#	--filterExpression 'ReadPosRankSum < -8.0' --filterName 'FAILS_HARD_FILTER_SNP_ReadPosRankSum' \
#	--variant {input.vcf} \
#	-o {output.vcf}
#	"""
	
############
rule vcf2tsv:
	input:
		OUTPUT_PATH+"/VCF/{sample}.vcf"
	output:
		OUTPUT_PATH+"/TSV/{sample}.tsv"
	params:
		script=WORKING_DIR+"/processVCF.R",
		generalRfun=WORKING_DIR+"/general.R"
	shell: """
		Rscript {params.script} {input} {output} {params.generalRfun}
	"""

rule mergeVCF:
	input:
		expand(OUTPUT_PATH+"/TSV/{sample}.tsv", sample=SAMPLES)
	output:
		OUTPUT_PATH+"/results/merged.tsv"
	shell: """
		Rscript {WORKING_DIR}/mergeVCF.R {OUTPUT_PATH} {WORKING_DIR}/general.R
	"""

rule MULTIQC:
	input:
		fastqc=expand(OUTPUT_PATH+"/fastQC/{sample}{pair}_fastqc{ext}", sample=SAMPLES,pair=PAIR_SUFFIX,ext=[".zip",".html"]),
		tsv=OUTPUT_PATH+"/results/merged.tsv"
	output: OUTPUT_PATH+"/results/multiqc_report.html"
	params:
		outpath = OUTPUT_PATH + "/results",
		cpu = 1
	shell: """
	multiqc -f -e general_stats -e tophat -e bowtie2 {OUTPUT_PATH} -o {params.outpath}
	"""
