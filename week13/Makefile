#
# Variant calling workflow.
#

# Accession number of the yeast genome.
ACC=GCF_000146045.2

# The reference file.
REF=refs/${ACC}.fa

# The GFF file.
GFF=refs/${ACC}.gff

# The sequencing read accession number.
SRR=SRR18902924

# The number of reads to get
N=5000

# The path to read 1
R1=reads/read1.fastq

# The path to read 2
R2=

# The resulting BAM file.
BAM=bam/HBR_1.bam

# The resulting variant VCF file (compressed!).
VCF=vcf/${SAMPLE}.vcf.gz

# Custom makefile settings.
SHELL = bash
.ONESHELL:
.SHELLFLAGS = -eu -o pipefail -c
.DELETE_ON_ERROR:
MAKEFLAGS += --warn-undefined-variables
MAKEFLAGS += --no-builtin-rules

# Print the usage of the makefile.
usage:
	@echo "#"
	@echo "# SNP call demonstration"
	@echo "#"
	@echo "# ACC=${ACC}"
	@echo "# SRR=${SRR}"
	@echo "# SAMPLE=${SAMPLE}"
	@echo "# BAM=${BAM}"
	@echo "# VCF=${VCF}"
	@echo "#"
	@echo "# make bam|vcf|all"
	@echo "#"

# Check that the bio toolbox is installed.
CHECK_FILE = src/run/genbank.mk
${CHECK_FILE}:
	@echo "#"
	@echo "# Please install toolbox with: bio code"
	@echo "#"
	@exit 1


# Index the genome
index:
# Get the reference genome and the annotations.
	make -f src/run/datasets.mk ACC=${ACC} REF=${REF} GFF=${GFF} run
# Download the sequence data.
	make -f src/run/sra.mk SRR=${SRR} R1=${R1} N=${N} run
# index
	make -f src/run/hisat2.mk index REF=${REF}

# Run the alignment for a single sample
alignment:
	make -f src/run/hisat2.mk \
			REF=${REF} R1=${R1} BAM=${BAM} \
			run

# Generate a bigwig coverage track from a BAM file
coverage:
	make -f src/run/wiggle.mk \
			REF=${REF} BAM=${BAM} \
			run

all: index alignment coverage 