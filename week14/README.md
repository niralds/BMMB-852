# Week 14: Perform a differential expression analysis

Your submission should include a readme, a makefile, and a design file.

You may re-use code and genomes from your previous submission.

Perform a differential expression analysis of a count matrix.

Take a count matrix, this count matrix may be one generated in a previous assignment or one that simulated.

Use method to identify genes/transcripts that show differential expression.

Draw a PCA plot and a heatmap for these genes.

Discuss the results. How many genes have you found. What kind of expression levels can you observe. How reliable does your data seem?

to generate the design.csv file 

````
bio search PRJNA831280 -H --csv > design.csv
````

this is followed by 

````
$ cat design.csv | head -10 | \
>     parallel --dry-run --lb -j 4 --colsep , --header : \
>     make all SRR={run_accession} SAMPLE={library_name}
````

the output is

````
make all SRR=SRR18902888 SAMPLE=GSM6062550
make all SRR=SRR18902889 SAMPLE=GSM6062549
make all SRR=SRR18902890 SAMPLE=GSM6062549
make all SRR=SRR18902891 SAMPLE=GSM6062548
make all SRR=SRR18902892 SAMPLE=GSM6062547
make all SRR=SRR18902894 SAMPLE=GSM6062546
make all SRR=SRR18902897 SAMPLE=GSM6062545
make all SRR=SRR18902901 SAMPLE=GSM6062542
make all SRR=SRR18902902 SAMPLE=GSM6062541
````

to run the makefile, type: 

````
make all 
````

the bigwig file looks like this: 

![alt text](<screenshot/Screenshot 2024-12-10 at 8.10.36 PM.png>)

it somewhat looks like the RNA seq worked

#### To process the count and design file, the makefile prodcues:
````
Initializing edgeR tibble dplyr tools ... done
# Tool: edgeR 
# Design: design.csv 
# Counts: counts.csv 
# Sample column: sample 
# Factor column: group 
# Factors: A B 
# Group A has 3 samples.
# Group B has 3 samples.
# Method: glm 
# Input: 20000 rows
# Removed: 15332 rows
# Fitted: 4668 rows
# Significant PVal:  413 ( 8.80 %)
# Significant FDRs:  137 ( 2.90 %)
# Results: edger.csv 
````

#### To evaluate the results, the following was the output 

````
# Tool: evaluate_results.r 
# 258 in counts.csv 
# 137 in edger.csv 
# 130 found in both
# 128 found only in counts.csv 
# 7 found only in edger.csv 
# Summary: summary.csv 
````

### The PCA plot 

![alt text](<plots/Screenshot 2024-12-15 at 9.06.29 PM.png>)

### The heatmap 

![alt text](<plots/Screenshot 2024-12-15 at 9.06.57 PM.png>)

#### The differentially expressed genes were

````
GENE-13123
GENE-6638
GENE-17107
GENE-4550
GENE-19918
GENE-8184
GENE-3463
GENE-3354
GENE-349
````

