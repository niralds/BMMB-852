# Set the trace
set -uex

# SRR number
SRR=SRR30805512

# Number of reads to sample
N=10000

# The output read names
R1=reads/${SRR}_1.fastq
R2=reads/${SRR}_2.fastq

# Trimmed read names
T1=reads/${SRR}_1.trimmed.fastq
T2=reads/${SRR}_2.trimmed.fastq

# The reads directory
RDIR=reads

# The reports directory
PDIR=reports

# ----- actions below ----

# Make the necessary directories
mkdir -p ${RDIR} ${PDIR}

# Download the FASTQ file
fastq-dump -X ${N} --split-files -O ${RDIR} ${SRR} 

# Run fastqc
fastqc -q -o ${PDIR} ${R1} ${R2}

# Run fastp and trim for quality
fastp --trim_poly_g \
      -i ${R1} -I ${R2} -o ${T1} -O ${T2} 

# Run fastqc
fastqc -q -o ${PDIR} ${T1} ${T2}

# Assumes that you have installed multiqc like so
# conda create -y -n menv multiqc 

# Run the multiqc from the menv environment on the reports directory
# The output of the tool is in the multiqc_report.html
micromamba run -n menv multiqc -o ${PDIR} ${PDIR}