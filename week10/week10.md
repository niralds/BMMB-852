# Week 10: Generate a variant call file

````
# Accession number
ACC = GCF_000146045.2

# The location of the reference file (FASTA format)
REF=refs/${ACC}
FAI=refs/${ACC}.fai

# SRR number
SRR=SRR11466177

# Index 
IDX_DIR = $(dir ${REF})idx
IDX = ${IDX_DIR}/$(notdir ${REF})
IDX_FILE = ${IDX}.ann

# Set the values for the read groups.
ID ?= run1
SM ?= sample1
LB ?= library1
PL ?= ILLUMINA

# Build the read groups tag.
RG = '@RG\tID:${ID}\tSM:${SM}\tLB:${LB}\tPL:${PL}'

# Number of reads to sample
N1=10000

# The output read names
R1=reads/${SRR}_1.fastq
R2=reads/${SRR}_2.fastq

# The resulting SAM/BAM files.
SAM=bam/${SRR}.sam
BAM=bam/${SRR}.bam
FTL=selectedqc.bam

# The bedgraph file.
BGR= ${BAM:.bam=.bedgraph}

# The wiggle file. Swap extension.
WIG= $(BAM:.bam=.bw)

# The variant file 
VCF ?= vcf/$(notdir $(basename ${BAM})).vcf.gz

# Additional bcf flags for pileup annotation.
PILE_FLAGS =  -d 100 --annotate 'INFO/AD,FORMAT/DP,FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/SP'

# Additional bcf flags for calling process.
CALL_FLAGS = --ploidy 2 --annotate 'FORMAT/GQ'



# Print the help
usage:
	@echo                "# ACC=${ACC}"
	@echo "make dir      # Make the directories "
	@echo "make ref        # Get the reference genome "
	@echo "make SRA              # Fetch sequence data from SRA (fastq format)"
	@echo "make index            # Index the genome"
	@ echo "make align           # Align the reads"
	@echo "make bam              # Sort the SAM file to BAM"
	@echo "make sam        # Index the BAM file"
	@echo "make unmap"     # reads that didn't align
	@echo "make stat"      # statistics of alignment
	@echo "make proper"    # properly-paired, on reverse strand, reads in first pair
	@echo "make qstat"     # statistics of filtered bam file


# Make the directories 
dir:
	mkdir -p bam reads refs

# Get the reference genome 
ref:
	bio fetch --format fasta ${ACC} > ${REF}

# Create the fasta index file.
refsam:
	samtools faidx ${REF}

# Fetch sequence data from SRA (fastq format)
SRA:
	fastq-dump -X ${N1} --outdir reads --split-files ${SRR} 

# Index the genome
index:
	mkdir ${IDX_DIR}
	bwa index -p ${IDX} ${REF}

# Align the reads
align:
	bwa mem -t 4 -R ${RG} ${IDX} ${R1} ${R2} > ${SAM}

# Sort the SAM file to BAM
bam file:
	cat ${SAM} | samtools sort > ${BAM}

# Index the BAM file
sam:
	samtools index ${BAM}

# Creates the bigwig file.
bigwig:
	# Create the directory.
	mkdir -p $(dir ${WIG})

	# Generate the temporary bedgraph file.
	LC_ALL=C; bedtools genomecov -ibam  ${BAM} -split -bg  | sort -k1,1 -k2,2n > ${BGR}
	
	# Convert the bedgraph file to bigwig.
	bedGraphToBigWig ${BGR} ${FAI} ${WIG}

# Generate the wiggle file.
run: ${WIG}
	@ls -lh ${WIG}

# Variant file 
VCF:
	mkdir vcf
	bcftools mpileup ${PILE_FLAGS} -O u -f ${REF} ${BAM} | \
		bcftools call ${CALL_FLAGS} -mv -O u | \
		bcftools norm -f ${REF} -d all -O u | \
		bcftools sort -O z > ${VCF}
````

Run the code 

````
make dir ref refsam SRA index align bam file sam bigwig run VCF 

````

Load the bam file and VCF file on IGV 

![alt text](<Screenshots/Screenshot 2024-11-10 at 7.29.51 PM.png>)

However, there is a variation show here, but the VCF file doesn't show there is a variation present 

![alt text](<Screenshots/Screenshot 2024-11-10 at 7.35.51 PM.png>)