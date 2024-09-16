# Reformat previous assignment 
[this is the link to the previous assignment
](https://github.com/niralds/BMMB-852/blob/main/week02/hwk_2.md)
# Visualize the GFF file - saccharomyces cerevisiae

#### To start, i activated conda by typing: 
>  conda activate bioinfo
#### Then i tyed the following commands: 
> datasets 

>  datasets summary genome

#### To get a summary of genome i picked - saccharomyces cerevisiae 

>  datasets summary genome accession GCF_000146045.2

This was found from the NCBI database 

#### The data was then passed through jq (jason quesry line)
>  datasets summary genome accession GCF_000146045.2 | jq

### Download the file 
>  datasets download genome accession GCF_000146045.2

#### unzip the file 
>  unzip ncbi_dataset.zip

this is what was in the file 
>  inflating: README.md               
  inflating: ncbi_dataset/data/assembly_data_report.jsonl  
  inflating: ncbi_dataset/data/GCF_000146045.2/GCF_000146045.2_R64_genomic.fna  
  inflating: ncbi_dataset/data/dataset_catalog.json  
  inflating: md5sum.txt

  #### this is what is inside the GCF_000146045.2_R64_genomic.fna file 
  >  cat ncbi_dataset/data/GCF_000146045.2/GCF_000146045.2_R64_genomic.fna | head
>NC_001133.9 Saccharomyces cerevisiae S288C chromosome I, complete sequence
ccacaccacacccacacacccacacaccacaccacacaccacaccacacccacacacacacatCCTAACACTACCCTAACACAGCCCTAATCTAACCCTGGCCAACCTGTCTCTCAACTTACCCTCCATTACCCTGCCTCCACTCGTTACCCTGTCCCATTCAACCATACCACTCCGAACCACCATCCATCCCTCTACTTACTACCACTCACCCACCGTTACCCTCCAATTACCCATATCCAACCCACTGCCACTTACCCTACCATTACCCTACCATCCACCATGACCTACTCACCATACTGTTCTTCTACCCACCATATTGAAACGCTAACAAATGATCGTAAATAACACACACGTGCTTACCCTACCACTTTATACCACCACCACATGCCATACTCACCCTCACTTGTATACTGATTTTACGTACGCACACGGATGCTACAGTATATACCATCTCAAACTTACCCTACTCTCAGATTCCACTTCACTCCATGGCCCATCTCTCACTGAATCAGTACCAAATGCACTCACATCATTATGCACGGCACTTGCCTCAGCGGTCTATACCCTGTGCCATTTACCCATAACGCCCATCATTATCCACATTTTGATATCTATATCTCATTCGGCGGTcccaaatattgtataaCTGCCCTTAATACATACGTTATACCACTTTTGCACCATATACTTACCACTCCATTTATATACACTTATGTC


#### we downloaded more exclusive data i.e genome, CDS, protein 
>   datasets download genome accession GCF_000146045.2 --include gff3,cds,protein,rna,genome

new files were added to the folder, that had to be unzipped
>.. unzip ncbi_dataset.zip 
Archive:  ncbi_dataset.zip
replace README.md? [y]es, [n]o, [A]ll, [N]one, [r]ename: A
  inflating: README.md               
  inflating: ncbi_dataset/data/assembly_data_report.jsonl  
  inflating: ncbi_dataset/data/GCF_000146045.2/GCF_000146045.2_R64_genomic.fna  
  inflating: ncbi_dataset/data/GCF_000146045.2/genomic.gff  
  inflating: ncbi_dataset/data/GCF_000146045.2/cds_from_genomic.fna  
  inflating: ncbi_dataset/data/GCF_000146045.2/protein.faa  
  inflating: ncbi_dataset/data/GCF_000146045.2/rna.fna  
  inflating: ncbi_dataset/data/dataset_catalog.json  
  inflating: md5sum.txt  

  #### genomic data visualization 

  download IGV for visualization 

  ## Separation of gff file into genes and coding sequence 
  #### To separate the genes
  >  cat ncbi_dataset/data/GCF_000146045.2/genomic.gff | awk ' $3=="gene" { print $0 } '

  the data shown is all genes as shown below 
  >  NC_001133.9	RefSeq	gene	1807	2169	.	-	.	ID=gene-YAL068C;Dbxref=GeneID:851229;Name=PAU8;end_range=2169,.;gbkey=Gene;gene=PAU8;gene_biotype=protein_coding;locus_tag=YAL068C;partial=true;start_range=.,1807
NC_001133.9	RefSeq	gene	2480	2707	.	+	.	ID=gene-YAL067W-A;Dbxref=GeneID:1466426;Name=YAL067W-A;end_range=2707,.;gbkey=Gene;gene_biotype=protein_coding;locus_tag=YAL067W-A;partial=true;start_range=.,2480
NC_001133.9	RefSeq	gene	7235	9016	.	-	.	ID=gene-YAL067C;Dbxref=GeneID:851230;Name=SEO1;end_range=9016,.;gbkey=Gene;gene=SEO1;gene_biotype=protein_coding;locus_tag=YAL067C;partial=true;start_range=.,7235

the data was added to its own file - gene.gff 
>  cat ncbi_dataset/data/GCF_000146045.2/genomic.gff | awk ' $3=="gene" { print $0 } ' > ncbi_dataset/data/GCF_000146045.2/gene.gff 

#### To separate the CDS 
>  cat ncbi_dataset/data/GCF_000146045.2/genomic.gff | awk ' $3=="CDS" { print $0 } ' > ncbi_dataset/data/GCF_000146045.2/CDS.gff

## creatinf a gff file 

Open visual studio code 
Paste the following on the second line 
>  chr1  .  .  1000  2000  .  .  .  .

There are spaces instead of tabs between each charecter. 
Remove the space manually and type in the search bar
>  ">indent using tabs" 

manually add the spaces by pressing tab after each character 

The chr1 will be the accession number of the file 

>  NC_001133.9	.	.	1000	2000	.	.	.	.

open the file on IGV and visualize the annotation 

#### creating strandedness (directionalitly)

duplicate the code on the second line, so that you have two annotations on IGV 

Make the changes as below 
>  NC_001133.9	.	.	1000	2000	.	+	.	.

>  NC_001133.9	.	.	1000	2000	.	-	.	.

There will be fish scales on the annotations 

Type the following:

NC_001133.9	.	CDS	1000	2000	.	+	.	Parent=transcript1;ID=cds1
NC_001133.9	.	CDS	2500	3500	.	-	.	Parent=transcript1;ID=cds2

This is what connects the two annotations 

