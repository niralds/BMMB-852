# Write a Script

In the previous lecture, you wrote a Markdown report. Rewrite all the code from that report as a Bash script. Your script should be:

Reusable: separate variable definitions from the code.
Documented: include comments that explain what each part of the script does.

Make a file 

````
code demo.sh 
````

The script is: 

````
# Set the error handling and trace
set -uex

# Define all variables at the top.

# The url is 
URL="https://ftp.ensembl.org/pub/current_gff3/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.112.gff3.gz"

# Name of gff3 file 
GFF="yeast.gff"

# ------ ALL THE ACTIONS FOLLOW ------

#Download gff3 file and store it in file 
curl ${URL} > ${GFF}.gz

#unzip the file 
gunzip ${GFF}.gz

#How many lines in this file 
cat ${GFF} | wc -l

# Make a new GFF file with only the features of type gene
cat ${GFF} | awk '$3 == "gene"' > genes.gff

# Print the number of genes
cat genes.gff | wc -l

# count how many top ten elements are in column 3 of gff3 file 
cat ${GFF} | cut -f 3 | sort | uniq -c | sort -rn | head 
````

Run the script 

````
bash demo.sh 
````

#### Reproducing someones data 

These are the variables i inputed: 

````
# Set the error handling and trace
set -uex

# Define all variables at the top.

# The url is 
URL="https://ftp.ensembl.org/pub/current_gff3/mus_musculus/Mus_musculus.GRCm39.112.gff3.gz"

# Name of gff3 file 
GFF="Mus_musculus.GRCm39.112"
````

My code worked and it gave an output 

I also added an additional code that wasn't in mine, but in theirs 

````
#clean up gff3 file
cat ${GFF} | grep -v '#' > clean.gff3
````

The code worked 


# Make use of ontologies 

Software installation 

````
# Update the bio package
pip install bio --upgrade

# Download the latest database
bio --download

# Use bio to explain a term
bio explain gene
````
Sequence Ontology

1. Choose a feature type from the GFF file and look up its definition in the sequence ontology.
2. Find both the parent terms and children nodes of the term.
3. Provide a short discussion of what you found.

Type the following: 

````
bio explain snorna
````

This is the output: 

A snoRNA (small nucleolar RNA) is any one of a class of
small RNAs that are associated with the eukaryotic nucleus
as components of small nucleolar ribonucleoproteins. They
participate in the processing or modifications of many RNAs,
mostly ribosomal RNAs (rRNAs) though snoRNAs are also known
to target other classes of RNA, including spliceosomal RNAs,
tRNAs, and mRNAs via a stretch of sequence that is
complementary to a sequence in the targeted RNA.

Parents:
- sncrna 
- snorna_primary_transcript (derives_from)

Children:
- c_d_box_snorna 
- h_aca_box_snorna 
- scarna 

Discussion 

SnoRNA help guide chemical modifications in other types of RNAs, so they could help stabilize other RNAs.