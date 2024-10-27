# Accession number
ACC = GCF_000146045.2

# The location of the reference file (FASTA format)
REF=refs/GCF_000146045.2-yeast.fa

# SRR number
SRR=SRR30805512

# Number of reads to sample
N1=10000

# The output read names
R1=reads/${SRR}_1.fastq
R2=reads/${SRR}_2.fastq

# The resulting SAM/BAM files.
SAM=bam/${SRR}.sam
BAM=bam/${SRR}.bam


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
	@echo "make bowtie           # using bowtie


# Make the directories 
dir:
	mkdir -p bam reads refs

# Get the reference genome 
ref:
	bio fetch --format fasta ${ACC} > ${REF}

# Fetch sequence data from SRA (fastq format)
SRA:
	fastq-dump -X ${N1} --outdir reads --split-files ${SRR} 

# Index the genome
index:
	bwa index ${REF}

# Align the reads
align:
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