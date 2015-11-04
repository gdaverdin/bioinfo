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

## Identifying gene models (gtf) that overlap with locations of differentially expressed genes (gtf)

I have placed a gtf file of differentially expressed exons called DifferentiallyExpressedExons.gtf in the directory. Use bedtools intersect to count the number of gene model exons that overlap with differentially expressed exons.

## Extracting alignments over specific intervals

You have competing collaborators who want your full dataset to study DPCK and a suite of other genes. You're not interested in sharing the entire dataset before publication, but still want to collaborate. You decide to share with them only the read alignments over the DPCK1 gene. How do you extract a region of aligned reads from a .bam alignment file?

Samtools can do this for you. The DPCK gene is start and stop codons are at 40716757 and 43112247, respectively. Let's isolate reads just for this region and name it DPCK1.sam. By default, samtools will output a SAM file. 

    samtools view Brca1Reads_aligned.sorted.dedup.realigned.bam chr17:40716757-43112247 > DPCK1.sam

    head DPCK1.sam

But SAM files are big, we should really compress this into .bam before saving it. What flag in samtools view can you add to output a .bam file instead of the default .sam file?

    samtools view -h
