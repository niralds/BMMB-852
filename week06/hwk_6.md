# Identify a bad sequencing dataset. You may need to evaluate multiple SRR numbers to find one with poor quality

I was looking for yeast sequecing data and these were the first few files i found, after analyzing the report 

SRR11466505_1.fastq

![alt text](<../../../Desktop/Screenshot 2024-10-05 at 1.17.43 PM.png>)

This data seems really, so i tried looking for more data, and found this file 
SRR30805512

![alt text](<../../../Desktop/Screenshot 2024-10-05 at 1.20.38 PM.png>)

So i went ahead and used this file to answer the following questions: 
Use data generated on the Illumina or Iontorrent platforms.
Write a script to download data from the SRA database.
````
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
fastp --cut_tail \
      -i ${R1} -I ${R2} -o ${T1} -O ${T2} 

# Run fastqc
fastqc -q -o ${PDIR} ${T1} ${T2}
````

Some of the red warnings went away, but there were still yellow warnings 

![alt text](<../../../Desktop/Screenshot 2024-10-05 at 1.24.43 PM.png>)

````
Run fastp and trim for quality
fastp --trim_poly_g \
      -i ${R1} -I ${R2} -o ${T1} -O ${T2} 
````

This is what i got: 

![alt text](<../../../Desktop/Screenshot 2024-10-05 at 1.28.10 PM.png>)

No change in the quality 

So i decided to just trim the poly G and this is what i got 

![alt text](<../../../Desktop/Screenshot 2024-10-05 at 1.30.18 PM.png>)

There are more green checkmarks and fewer yellow warnings 