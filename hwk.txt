1. Tell us a bit about the organism
I chose saccharomyces_cerevisiae because I want to learn more about them. According to ChatGpt : Saccharomyces cerevisiae is a species of yeast that is widely used in baking, brewing, and scientific research. It's a single-celled fungus that's well-known for its role in fermenting sugars to produce carbon dioxide and alcohol. This process is crucial in the production of bread, beer, and wine. In the baking process, the carbon dioxide produced by the yeast makes dough rise, while in brewing and winemaking, it contributes to the alcohol content and flavor.

Beyond its practical applications, Saccharomyces cerevisiae is also an important model organism in biological research. It's used to study cellular processes, genetics, and biochemistry due to its relatively simple eukaryotic cell structure and well-understood genetics. Its genetic and molecular pathways are often conserved in other eukaryotes, including humans, making it a valuable tool for scientific discovery.

2. How many features does the file contain?
there are 8 features:
Line 1: I	R64-1-1	chromosome	1	230218	.	.	.	ID=chromosome:I;Alias=BK006935.2
Line 2: I	sgd	gene	335	649	.	+	.	ID=gene:YAL069W;biotype=protein_coding;description=Dubious open reading frame%3B unlikely to encode a functional protein%2C based on available experimental and comparative sequence data [Source:SGD%3BAcc:S000002143];gene_id=YAL069W;logic_name=sgd
(bioinfo) 

3. How many sequence regions (chromosomes) does the file contain?
there are 17 chromosomes 
nshah@Nirals-MacBook-Air ~/work/edu/unixtut
$ cat clean.gff | cut -f 1 | sort | uniq
I
II
III
IV
IX
Mito
V
VI
VII
VIII
X
XI
XII
XIII
XIV
XV
XVI
(bioinfo) 
nshah@Nirals-MacBook-Air ~/work/edu/unixtut
$ cat clean.gff | cut -f 1 | sort | uniq | wc
      17      17      59
(bioinfo) 

4. How many genes are listed for this organism? 
7125 genes 
nshah@Nirals-MacBook-Air ~/work/edu/unixtut
$ cat clean.gff | cut -f 3 | grep gene | wc -l
    7125
(bioinfo) 

5. What are the top-ten most annotated feature types (column 3) across the genome?
nshah@Nirals-MacBook-Air ~/work/edu/unixtut
$ cat clean.gff | cut -f 3 | sort | uniq -c | sort -rn | head
7507 exon
6913 CDS
6600 mRNA
6598 gene
 424 ncRNA_gene
 299 tRNA
  91 transposable_element_gene
  91 transposable_element
  77 snoRNA
  24 rRNA

6. Having analyzed this GFF file, does it seem like a complete and well-annotated organism?
yeah it seems to be well annotated 

7. Share any other insights you might note.
I kept the downloaded file open while completing this assignment to really understand what we were doing and where we were pulling information from 
