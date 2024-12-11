# Week 13: Generate an RNA-Seq count matrix

Your submission should include a readme, a makefile, and a design file.

You may re-use code and genomes from your previous submission. You will need a genome and a GTF/GFF annotation reference. You may also perform a transcriptome-based quantification if you choose.

For your genome of interest, find RNA-Seq data at SRA. Most RNA-Seq projects include multiple datasets.

Select at least three SRR datasets from the same project.

Write code to perform an RNA-Seq analysis on each dataset. The final result of your code should be a count matrix that summarizes read counts for each dataset.

Include screenshots from IGV that prove that your data is indeed RNA-Seq data. 

Discuss a few lines of the resulting count matrix. Visually identify rows where the counts show consistent gene expression levels.

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
cat design.csv | head -25 | \
    parallel --lb -j 4 --colsep , --header : \
    make all SRR={run_accession} SAMPLE={library_name}

````

however, when i run this on my makefile, it gives me an error 

If i run the SRR number individually, it works and makes a bigwig file. i.e 

````
make all 
````

the bigwig file looks like this: 

![alt text](<screenshot/Screenshot 2024-12-10 at 8.10.36 PM.png>)

it somewhat looks like the RNA seq worked