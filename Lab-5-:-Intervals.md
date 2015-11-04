Today we will go over a set of practical exercises regarding intervals in genomic data. Intervals are critical for genomics; every gene annotation, transposon annotation, and read alignment is really just an interval of locations along a genome or transcriptome -- a start and stop location. 

Today's lab will focus largely on using the bedtools and samtools suites of programs. We will use the data and your outputs from the Variant Calling lab. I have placed these files here so you can copy them into your own directory:
    
    /home/student/binf4550/data/05.Intervals

Let's see what bedtools can do:

    /usr/local/bedtools/latest/bin/bedtools -h

## Identify gene models (gtf) that intersect with variant calls (vcf)

Goal: Which gene models intersect with a SNP that we called?

We'll start by using bedtools intersect to identify the intersection of two or more annotations.

![](http://bedtools.readthedocs.org/en/latest/_images/intersect-glyph.png)
Source: http://quinlanlab.org/tutorials/cshl2013/bedtools.html

    /usr/local/bedtools/latest/bin/bedtools intersect -h
    /usr/local/bedtools/latest/bin/bedtools intersect -a Hg19.Chr17.UCSC-3.exons.gtf -b Brca_raw_variants.vcf

And look -- there is the answer to last week's lab homework -- the only SNP located in an exon. Based on the above picture, add the flag "-wa" and see what happens. Then "-wb" instead. 

If you wanted to identify every gene model that DIDNT overlap with a SNP, then add "-v". 