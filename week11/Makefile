# Accession number
ACC = GCF_000146045.2

# The number of reads
N=600000

GENOMIC=${ACC}.fa

# The location of the reference file (FASTA format)
REF=${ACC}.fa
FAI=${ACC}.fa.fai

# SRR number
SRR=SRR11467335

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

# Lengh of the reads
L=100

# The reads directory
RDIR=reads

# The reports directory
PDIR=reports

# The files to write the reads to
DS=reads
R1=reads/wgsim_read1.fq
R2=reads/wgsim_read2.fq

# The output read names
R11=reads/${SRR}_1.fastq
R12=reads/${SRR}_2.fastq

# Trimmed read names
T1=reads/${SRR}_1.trimmed.fastq
T2=reads/${SRR}_2.trimmed.fastq

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


# download the genome file 
genome: 
	# Download the genome
	datasets download genome accession ${ACC} --include genome,gff3,gtf,rna,protein,cds
	# Unpack the data. (Overwrite files if they already exist)
	unzip -o ncbi_dataset.zip
	# Make a link to a simpler name
	ln -sf ncbi_dataset/data/GCF_000146045.2/GCF_000146045.2_R64_genomic.fna ${ACC}.fa
	@echo genome downloaded


#simulate reads for the genome
simulate:
	# Make the directory that will hold the reads extracts 
	# the directory portion of the file path from the read
	mkdir -p ${DS}
	# Simulate with no errors and no mutations
	wgsim -N ${N} -1 ${L} -2 ${L} -r 0 -R 0 -X 0 ${GENOMIC} ${R1} ${R2}
	# Run read statistics
	seqkit stats ${R1} ${R2}
	@echo reads stimulated 

# dowmload reads from SRA
download: 
	# Make the necessary directories
	mkdir -p ${RDIR} ${PDIR}

	# Download the FASTQ file
	fastq-dump -X ${N1} --outdir reads --split-files ${SRR} 
	@echo reads dowmloaded from SRA

# trim reads
report: 
	# Run fastqc
	fastqc -q -o ${PDIR} ${R11} ${R12}
trim:
	# Run fastp and trim for quality
	fastp --trim_poly_g \
		-i ${R11} -I ${R12} -o ${T1} -O ${T2} 

	# Run fastqc
	fastqc -q -o ${PDIR} ${T1} ${T2}
	@echo reads trimmed

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
	mkdir bam
	bwa mem -t 4 ${IDX} ${R1} ${R2} > ${SAM}

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

# VEP needs a sorted and compressed GFF file.
${GFF}.gz:
	# Sort and compress the GFF file
	# Needs the double $ to pass the $ from make to bash
	cat ${GFF} | sort -k1,1 -k4,4n -k5,5n -t$$'\t' | bgzip -c > ${GFF}.gz

	# Index the GFF file
	tabix -p gff ${GFF}.gz

# VEP is installed in the environment called vep
vep: ${GFF}.gz
	mkdir -p results
	micromamba run -n vep \
		~/src/ensembl-vep/vep \
		-i ${VCF} \
		-o results/vep.txt \
		--gff ${GFF}.gz \
		--fasta ${REF} \
		--force_overwrite 

	# Show the resulting files
	ls -lh results/*

# Cleanup the downloaded files.
clean:
	rm -rf ncbi_dataset/data/${ACC}
	rm -f md5sum.txt ncbi_dataset.zip
	rm -f fastp.html fastp.json
	rm -f README.md

	# Mark the targets that do not create files.
.PHONY: usage genome download clean