# Write a makefile 

This markdown reports is about how to simulate reads obtain and trim reads for a realistic dataset using a makefile. 

First create a file called:

````
Makefile
````

This is the file: 

````
# Accession number
ACC = GCF_000146045.2

# The location of the genome file (FASTA format)
GENOME=yeast.fa

# The number of reads
N=600000

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

# Print the help
usage:
	@echo                "# ACC=${ACC}"
	@echo "make info      # summary information on the genome"
	@echo "make genome    # download the genome file"
	@echo "make stimulate # stimulate reads for the genome"
	@echo "make donwload  # download reads from SRA" 
	@ echo "make trim     # trim reads"
	@echo "make clean     # remove the downloaded files"

# Print the summary information on the genome.
info:
	@datasets summary genome accession ${ACC} | jq
	@echo info provided

# download the genome file 
genome: 
	# Download the genome
	datasets download genome accession ${ACC}
	# Unpack the data. (Overwrite files if they already exist)
	unzip -o ncbi_dataset.zip
	# Make a link to a simpler name
	ln -sf ncbi_dataset/data/GCF_000146045.2/GCF_000146045.2_R64_genomic.fna yeast.fa
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
	fastq-dump -X ${N1} --split-files -O ${RDIR} ${SRR} 
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
````

Different parts of the file can be run separately 

````
make usage
````

Outputs: 

````
# ACC=GCF_000146045.2
make info      # summary information on the genome
make genome    # download the genome file
make stimulate # stimulate reads for the genome
make donwload  # download reads from SRA
make trim     # trim reads
make clean     # remove the downloaded files
````

To run everything at once: 

````
make all 
````

to clean the directory: 

````
make clean
````