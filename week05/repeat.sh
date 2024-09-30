# The location of the genome file (FASTA format)
GENOME=yeast.fa

# The number of reads
N=3000

# Lengh of the reads
L=200

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