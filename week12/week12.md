# Week 12: Automate a VCF calling pipeline 

Create a Makefile that can produce a VCF file by downloading a reference genome, indexing it, and then downloading fastq files from SRA, aligning them to the reference and calling variants.
Create a README.md file that explains how to run the Makefile
Collect a set of samples from the SRA database that match your genome.
Create a design.csv file that lists the samples to be processed.
Using GNU parallel or any other method of your choice run the Makefile on all (or a a subset) of the samples
Merge the resulting VCF files into a single one.
Discuss what you see in your VCF file.

to generate the design.csv file 

````
bio search PRJNA956066 -H --csv > design.csv
````

this is followed by 

````
$ cat design.csv | head -10 | \
>     parallel --dry-run --lb -j 4 --colsep , --header : \
>     make all SRR={run_accession} SAMPLE={library_name}
````

the output is

````
make all SRR=SRR24181940 SAMPLE=1409192300
make all SRR=SRR24181942 SAMPLE='414914105;413976152'
make all SRR=SRR24181943 SAMPLE='488145053;487195227'
make all SRR=SRR24181945 SAMPLE='576049962;574592921'
make all SRR=SRR24181948 SAMPLE=3088340800
make all SRR=SRR24181951 SAMPLE=3069795700
make all SRR=SRR24181962 SAMPLE='637626739;636884693'
make all SRR=SRR24181965 SAMPLE=2123144500
make all SRR=SRR24181966 SAMPLE='584909790;584031030'
`````

to run the makefile and create all the merged VCF file: 

````
cat design.csv | head -25 | \
    parallel --lb -j 4 --colsep , --header : \
    make all SRR={run_accession} SAMPLE={library_name}

````

this is what the merged VCF files look like 

![alt text](<images/Screenshot 2024-12-03 at 3.17.23 PM.png>)
