# Accession number
ACC = GCF_000146045.2

# The number of reads
N=600000

GENOME=${ACC}.fa

# Lengh of the reads
L=100

# The files to write the reads to
DS=reads
R1=reads/wgsim_read1.fq
R2=reads/wgsim_read2.fq

# SRR number
SRR=SRR30805512

# Number of reads to sample
N1=10000

# The output read names
R11=reads/${SRR}_1.fastq
R12=reads/${SRR}_2.fastq

# Trimmed read names
T1=reads/${SRR}_1.trimmed.fastq
T2=reads/${SRR}_2.trimmed.fastq

# The reads directory
RDIR=reads

# The reports directory
PDIR=reports

# The resulting SAM/BAM files.
REF=${ACC}.fa
SAM=bam/${SRR}.sam
BAM=bam/${SRR}.bam

# Print the help
usage:
	@echo                "# ACC=${ACC}"
	@echo "make info      # summary information on the genome"
	@echo "make genome    # download the genome file"
	@echo "make stimulate # stimulate reads for the genome"
	@echo "make donwload  # download reads from SRA" 
	@ echo "make trim     # trim reads"
	@ echo "make SRA      # fetch sequence data from SRA"
	@echo "make index            # Index the genome"
	@ echo "make align           # Align the reads"
	@echo "make bam file              # Sort the SAM file to BAM"
	@echo "make sam        # Index the BAM file"
	@echo "make stats      # generate alignment statistics"
	@echo "make clean     # remove the downloaded files"

# Print the summary information on the genome.
info:
	@datasets summary genome accession ${ACC} | jq
	@echo info provided

# download the genome file 
genome: 
	# Download the genome
	mkdir refs
	datasets download genome accession ${ACC}
	# Unpack the data. (Overwrite files if they already exist)
	unzip -o ncbi_dataset.zip
	# Make a link to a simpler name
	ln -sf ncbi_dataset/data/GCF_000146045.2/GCF_000146045.2_R64_genomic.fna ${ACC}.fa
	@echo genome downloaded


#stimulate reads for the genome
stimulate:
	# Make the directory that will hold the reads extracts 
	# the directory portion of the file path from the read
	mkdir -p ${DS}
	# Simulate with no errors and no mutations
	wgsim -N ${N} -1 ${L} -2 ${L} -r 0 -R 0 -X 0 ${GENOME} ${R1} ${R2}
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
trim: 
	# Run fastqc
	fastqc -q -o ${PDIR} ${R11} ${R12}

	# Run fastp and trim for quality
	fastp --trim_poly_g \
		-i ${R11} -I ${R12} -o ${T1} -O ${T2} 

	# Run fastqc
	fastqc -q -o ${PDIR} ${T1} ${T2}
	@echo reads trimmed

# Fetch sequence data from SRA (fastq format)
SRA:
	fastq-dump -X ${N1} --outdir reads --split-files ${SRR} 

# Index the genome
index:
	bwa index ${REF}

# Align the reads
align:
	mkdir bam
	bwa mem -t 4 ${REF} ${R1} ${R2} > ${SAM}

# Sort the SAM file to BAM
bam file:
	cat ${SAM} | samtools sort > ${BAM}

# Index the BAM file
sam:
	samtools index ${BAM}

# Generate alignment statistics
stats:
	samtools flagstat ${BAM}
	samtools flagstat ${R1}
	samtools flagstat ${R2}

# Cleanup the downloaded files.
clean:
	rm -rf ncbi_dataset/data/${ACC}
	rm -f md5sum.txt ncbi_dataset.zip
	rm -f fastp.html fastp.json
	rm -f README.md

	# Mark the targets that do not create files.
	.PHONY: usage genome download clean