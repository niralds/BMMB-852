# Set the error handling and trace
set -uex

# Define all variables at the top.

# The url is 
URL="https://ftp.ensembl.org/pub/current_gff3/mus_musculus/Mus_musculus.GRCm39.112.gff3.gz"

# Name of gff3 file 
GFF="Mus_musculus.GRCm39.112"

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

#clean up gff3 file
cat ${GFF} | grep -v '#' > clean.gff3

