#
# Variant calling workflow.
#

# Accession number of the yeast genome.
ACC=GCF_000146045.2

# The reference file.
REF=refs/${ACC}.fa

# The GFF file.
GFF=refs/${ACC}.gff

# The GTF file
GTF=refs/${ACC}.gtf

# The counts in tab delimited format.
COUNTS_TXT = res/counts.txt

# Final combinted counts in CSV format.
COUNTS = res/counts.csv

# The design file
DESIGN = design.csv


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

# Generate count matrix 
count: 
	# Make the directory name for the counts
	mkdir res

	# Count the features
	featureCounts -a ${GFF} -o ${COUNTS_TXT} ${BAM}

# The final counts in CSV format.
	micromamba run -n stats Rscript src/r/format_featurecounts.r -c ${COUNTS_TXT} -o ${COUNTS}

# Trigger the counting explicitly
	@ls -lh ${COUNTS_TXT}
	@ls -lh ${COUNTS}

simulate: 
	micromamba run -n stats Rscript src/r/simulate_counts.r
	micromamba run -n stats Rscript src/r/edger.r

evaluate: 
	micromamba run -n stats Rscript  src/r/evaluate_results.r  -a counts.csv -b edger.csv

	cat edger.csv | cut -f 1 -d ,  | head -10

	micromamba run -n stats src/r/plot_pca.r -c edger.csv 
	micromamba run -n stats src/r/plot_heatmap.r -c edger.csv

all: index alignment coverage count simulate evaluate