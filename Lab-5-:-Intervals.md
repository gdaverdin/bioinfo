Today we will go over a set of practical exercises regarding intervals in genomic data. Intervals are critical for genomics; every gene annotation, transposon annotation, and read alignment is really just an interval of locations along a genome or transcriptome -- a start and stop location. 

Today's lab will focus largely on using the bedtools suite of programs. We will use the data and your outputs from the Variant Calling lab. I have placed these files here so you can copy them into your own directory:
    
    /home/student/binf4550/data/05.Intervals

Let's see what bedtools can do:

    /usr/local/bedtools/latest/bin/bedtools -h

## Identify gene models (gtf) that intersect with variant calls (vcf)

Goal: Which gene models intersect with a SNP that we called?

    /usr/local/bedtools/latest/bin/bedtools intersect -h
    /usr/local/bedtools/latest/bin/bedtools intersect -a Hg19.Chr17.UCSC-3.gtf -b Brca_raw_variants.vcf

    /usr/local/bedtools/latest/bin/bedtools intersect -a Hg19.Chr17.UCSC-3.gtf -b Brca_raw_variants.vcf
    chr17	hg19_spAnnot	CDS	43107539	43107539	1000.000000	.	2	gene_id "DPCK"; transcript_id "DPCK_dup1";
    chr17	hg19_spAnnot	exon	43107539	43107539	1000.000000	.	.	gene_id "DPCK"; transcript_id "DPCK_dup1";

And look -- there is the answer to last week's lab homework -- the only SNP located in an exon.