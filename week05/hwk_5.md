# Sequencing Technologies

Select a genome, then download the corresponding FASTA file.

````
datasets download genome accession GCF_000146045.2
````

What is inside the file? 

````
cat GCF_000146045.2_R64_genomic.fna | head
````

This is the output it gives: 

````
>NC_001133.9 Saccharomyces cerevisiae S288C chromosome I, complete sequence
ccacaccacacccacacacccacacaccacaccacacaccacaccacacccacacacacacatCCTAACACTACCCTAAC
ACAGCCCTAATCTAACCCTGGCCAACCTGTCTCTCAACTTACCCTCCATTACCCTGCCTCCACTCGTTACCCTGTCCCAT
TCAACCATACCACTCCGAACCACCATCCATCCCTCTACTTACTACCACTCACCCACCGTTACCCTCCAATTACCCATATC
CAACCCACTGCCACTTACCCTACCATTACCCTACCATCCACCATGACCTACTCACCATACTGTTCTTCTACCCACCATAT
TGAAACGCTAACAAATGATCGTAAATAACACACACGTGCTTACCCTACCACTTTATACCACCACCACATGCCATACTCAC
CCTCACTTGTATACTGATTTTACGTACGCACACGGATGCTACAGTATATACCATCTCAAACTTACCCTACTCTCAGATTC
CACTTCACTCCATGGCCCATCTCTCACTGAATCAGTACCAAATGCACTCACATCATTATGCACGGCACTTGCCTCAGCGG
TCTATACCCTGTGCCATTTACCCATAACGCCCATCATTATCCACATTTTGATATCTATATCTCATTCGGCGGTcccaaat
attgtataaCTGCCCTTAATACATACGTTATACCACTTTTGCACCATATACTTACCACTCCATTTATATACACTTATGTC
````

#### Size of the file 

````
ls -lh GCF_000146045.2_R64_genomic.fna 
````

Output: 

````
12M Sep 15 16:28 GCF_000146045.2_R64_genomic.fna
````

So it is 12MB 

#### The total size of the genome

````
datasets summary genome accession GCF_000146045.2 | jq
````

After reading through the output, i saw this:

````
"total_sequence_length": "12071326"
````
This makes sense because budding yeast is 12 million in size. 

#### The number of chromosomes in the genome

From the same code i ran above, i found: 

````
"total_number_of_chromosomes": 16,
````

#### The name(id) and length of each chromosome in the genome

````
awk '/^>/ { if (seqlen) {
>               print seqlen
>               }
>             print
> 
>             seqtotal+=seqlen
>             seqlen=0
>             seq+=1
>             next
>             }
>     {
>     seqlen += length($0)
>     }     
>     END{print seqlen
>         print seq" sequences, total length " seqtotal+seqlen
>     }' yeast.fa

````

Output: 

````
>NC_001133.9 Saccharomyces cerevisiae S288C chromosome I, complete sequence
230218
>NC_001134.8 Saccharomyces cerevisiae S288C chromosome II, complete sequence
813184
>NC_001135.5 Saccharomyces cerevisiae S288C chromosome III, complete sequence
316620
>NC_001136.10 Saccharomyces cerevisiae S288C chromosome IV, complete sequence
1531933
>NC_001137.3 Saccharomyces cerevisiae S288C chromosome V, complete sequence
576874
>NC_001138.5 Saccharomyces cerevisiae S288C chromosome VI, complete sequence
270161
>NC_001139.9 Saccharomyces cerevisiae S288C chromosome VII, complete sequence
1090940
>NC_001140.6 Saccharomyces cerevisiae S288C chromosome VIII, complete sequence
562643
>NC_001141.2 Saccharomyces cerevisiae S288C chromosome IX, complete sequence
439888
>NC_001142.9 Saccharomyces cerevisiae S288C chromosome X, complete sequence
745751
>NC_001143.9 Saccharomyces cerevisiae S288C chromosome XI, complete sequence
666816
>NC_001144.5 Saccharomyces cerevisiae S288C chromosome XII, complete sequence
1078177
>NC_001145.3 Saccharomyces cerevisiae S288C chromosome XIII, complete sequence
924431
>NC_001146.8 Saccharomyces cerevisiae S288C chromosome XIV, complete sequence
784333
>NC_001147.6 Saccharomyces cerevisiae S288C chromosome XV, complete sequence
1091291
>NC_001148.4 Saccharomyces cerevisiae S288C chromosome XVI, complete sequence
948066
>NC_001224.1 Saccharomyces cerevisiae S288c mitochondrion, complete genome
85779
17 sequences, total length 12157105
````


# Generate a simulated FASTQ output for a sequencing instrument of your choice.  Set the parameters so that your target coverage is 10x.

#### Generate reads using wgsim 

yeast genome = 12 million 
number of reads = 100 
10 X coverage for paired end means 600,000 reads

Code: 

````
#
# Simulate reads with wgsim
#

# Set the trace
set -uex

# The location of the genome file (FASTA format)
GENOME=GCF_000146045.2_R64_genomic.fna

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
````

To run code: 

````
bash hwk_5.sh
````
Output: 

````
New version of client (16.30.1) available at https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/LATEST/mac/datasets.
Collecting 1 genome record [================================================] 100% 1/1
Downloading: ncbi_dataset.zip    3.85MB valid zip structure -- files not checked
Validating package [================================================] 100% 5/5
Archive:  ncbi_dataset.zip
  inflating: README.md               
  inflating: ncbi_dataset/data/assembly_data_report.jsonl  
  inflating: ncbi_dataset/data/GCF_000146045.2/GCF_000146045.2_R64_genomic.fna  
  inflating: ncbi_dataset/data/dataset_catalog.json  
  inflating: md5sum.txt              
+ GENOME=yeast.fa
+ N=600000
+ L=100
+ R1=reads/wgsim_read1.fq
+ R2=reads/wgsim_read2.fq
++ dirname reads/wgsim_read1.fq
+ mkdir -p reads
+ wgsim -N 600000 -1 100 -2 100 -r 0 -R 0 -X 0 yeast.fa reads/wgsim_read1.fq reads/wgsim_read2.fq
[wgsim] seed = 1727724825
[wgsim_core] calculating the total length of the reference sequence...
[wgsim_core] 17 sequences, total length: 12157105
+ seqkit stats reads/wgsim_read1.fq reads/wgsim_read2.fq
processed files:  2 / 2 [======================================] ETA: 0s. done
file                  format  type  num_seqs     sum_len  min_len  avg_len  max_len
reads/wgsim_read1.fq  FASTQ   DNA    600,000  60,000,000      100      100      100
reads/wgsim_read2.fq  FASTQ   DNA    600,000  60,000,000      100      100      100
````

How many reads have you generated?
600,000 sequences gives a total of 60,000,000 read length, which is 120,000,000 from paired end. 

What is the average read length?
100 

How big are the FASTQ files? 143 Mb 

Compress the file 

````
gzip wgsim_read1.fq 
````

File size: 28Mb 

#### Discuss whether you could get the same coverage with different parameter settings (read length vs. read number).

Changing the parameters i.e The number of reads =300000, and lengh of the reads =200 gives ~ 60,000,000 reads 

So, same coverage is found 

# How much data would be generated when covering the Yeast,  the Drosophila or the Human genome at 30x?


yeast genome = 12 MB
number of reads = 100 
paireed end, 30 X coverage for paired end means 1,800,000 reads
Each FastQ before compression: 360MB 
Each astQ file after compression: 72MB

Human genome = 3.1 GB 
number of reads = 100 
paired end, 30 X coverage for paired end means 465,000,000 reads
Each FastQ before compression: 93GB
Each FastQ file after compression: 18.6GB 

Drosophila genome = 143 MB
number of reads = 100 
paireed end, 30 X coverage for paired end means 21,450,000 reads
Each FastQ before compression: 4.29GB 
Each astQ file after compression: 858MB




