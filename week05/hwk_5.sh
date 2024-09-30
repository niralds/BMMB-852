# Download the genome
datasets download genome accession GCF_000146045.2

# Unpack the data. (Overwrite files if they already exist)
unzip -o ncbi_dataset.zip

# Make a link to a simpler name
ln -sf ncbi_dataset/data/GCF_000146045.2/GCF_000146045.2_R64_genomic.fna yeast.fa

#
# Simulate reads with wgsim
#

# Set the trace
set -uex

# The location of the genome file (FASTA format)
GENOME=yeast.fa

# The number of reads
N=600000

# Lengh of the reads
L=100

# The files to write the reads to
R1=reads/wgsim_read1.fq
R2=reads/wgsim_read2.fq

# --- Simulation actions below ---

# Make the directory that will hold the reads extracts 
# the directory portion of the file path from the read
mkdir -p $(dirname ${R1})

# Simulate with no errors and no mutations
wgsim -N ${N} -1 ${L} -2 ${L} -r 0 -R 0 -X 0 \
      ${GENOME} ${R1} ${R2}

# Run read statistics
seqkit stats ${R1} ${R2}