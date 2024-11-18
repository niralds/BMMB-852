# Week 11: Establish the effects of variants

This assignment requires writing a Makefile and a markdown report.

You can reuse the Makefile developed for your previous assignment, which generated a VCF file.

Evaluate the effects of the variants in your VCF file.

Try using a software tool like VEP or snpEff.  Add the effect prediction steps to your Makefile and make them part of the workflow.

If, for some reason, you can't make any of the variant effect prediction software work, you may use visual inspection via IGV to describe the effect of variants.

Find variants with different effects.

Write a markdown report that summarizes the process and your results.

#### run the code in the makefile 

When i type

````
make vep
````

It gives the following output 

````
$ make vep
# Sort and compress the GFF file
# Needs the double to pass the from make to bash
cat  | sort -k1,1 -k4,4n -k5,5n -t$'\t' | bgzip -c > .gz
````` 

doesn't load any further than this 

So i defualted to using the website 

This is the screenshot from IGV that shows the bam file, GFF file and VCF file. Notice the variant, and that there is a C mutation and the CDS track shows an A. And in the second screenshot, when i click the variant, it shows there is a A to C mutation 

![alt text](<screenshots/Screenshot 2024-11-17 at 8.52.05 PM.png>)

![alt text](<screenshots/Screenshot 2024-11-17 at 8.52.15 PM.png>)


This is the screenshot from the website 

![alt text](<screenshots/Screenshot 2024-11-17 at 8.46.34 PM.png>)