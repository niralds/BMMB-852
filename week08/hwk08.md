# Generate a BAM alignment file

First create a file called:

````
Makefile
````

use this code in the makefile:
````
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
	@echo "make refs        # Get the reference genome "
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
refs:
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
````

To run everything at once: 

````
make all 
````

in the bam directory, there should be a two bam files

Open these in IGV

#### Visualize the resulting BAM files for your simulated reads and the reads downloaded from SRA.

![alt text](<../../BMMB-852/week08/Screenshot 2024-10-20 at 5.47.24 PM.png>)

#### Generate alignment statistics for the reads from both sources, simulated and SRA.

this is the output in the terminal 

````
samtools flagstat bam/SRR30805512.bam
20108 + 0 in total (QC-passed reads + QC-failed reads)
20000 + 0 primary
0 + 0 secondary
108 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
363 + 0 mapped (1.81% : N/A)
255 + 0 primary mapped (1.27% : N/A)
20000 + 0 paired in sequencing
10000 + 0 read1
10000 + 0 read2
0 + 0 properly paired (0.00% : N/A)
20 + 0 with itself and mate mapped
235 + 0 singletons (1.17% : N/A)
8 + 0 with mate mapped to a different chr
3 + 0 with mate mapped to a different chr (mapQ>=5)


samtools flagstat reads/SRR30805512_1.fastq
10000 + 0 in total (QC-passed reads + QC-failed reads)
10000 + 0 primary
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
0 + 0 mapped (0.00% : N/A)
0 + 0 primary mapped (0.00% : N/A)
0 + 0 paired in sequencing
0 + 0 read1
0 + 0 read2
0 + 0 properly paired (N/A : N/A)
0 + 0 with itself and mate mapped
0 + 0 singletons (N/A : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)


samtools flagstat reads/SRR30805512_2.fastq
10000 + 0 in total (QC-passed reads + QC-failed reads)
10000 + 0 primary
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
0 + 0 mapped (0.00% : N/A)
0 + 0 primary mapped (0.00% : N/A)
0 + 0 paired in sequencing
0 + 0 read1
0 + 0 read2
0 + 0 properly paired (N/A : N/A)
0 + 0 with itself and mate mapped
0 + 0 singletons (N/A : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)
````