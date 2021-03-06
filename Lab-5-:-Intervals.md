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

## Find the closest SNP to each exon

I would like to know -- for each gene model, what is the closest SNP that we called?

    /usr/local/bedtools/latest/bin/bedtools closest -a Hg19.Chr17.UCSC-3.gtf -b Brca_raw_variants.vcf

Weird. It's throwing an error saying that our .gtf isn't sorted. 

    Error: Sorted input specified, but the file Hg19.Chr17.UCSC-3.exons.gtf has the following out of order record
    chr17	hg19_spAnnot	exon	131559	131645	1000.000000	.	.	gene_id "zinc finger"; transcript_id "zinc finger";

OK. bedtools can sort .gtf files for us. Here's what I did to fix this:

    /usr/local/bedtools/latest/bin/bedtools sort -i Hg19.Chr17.UCSC-3.exons.gtf > Hg19.Chr17.UCSC-3.exons.sort.gtf

Then we can run bedtools closest again. 

    /usr/local/bedtools/latest/bin/bedtools closest -a Hg19.Chr17.UCSC-3.exons.sort.gtf -b Brca_raw_variants.vcf

Redirect the output using > to a file of your choosing. I would probably name it something like "Hg19.Brca_variants.closest.out"

## Identifying gene models (gtf) that overlap with locations of differentially expressed genes (gtf)

I have placed a gtf file of differentially expressed exons called DifferentiallyExpressedExons.gtf in the directory. Use bedtools intersect to count the number of gene model exons that overlap with differentially expressed exons.

## Calculating the read coverage at specific genes

Since this dataset was a gene capture enrichment, we might expect variable coverage at multiple genes in the genome. I would like to know the average read coverage across DPCK. 

    /usr/local/bedtools/latest/bin/bedtools coverage -a Hg19.Chr17.UCSC-3.exons.gtf -b Brca1Reads_aligned.sorted.dedup.realigned.bam | grep "DPCK"

(See how I piped the output into grep to ONLY select exons from the DPCK gene?)

    Default Output:
	 After each entry in A, reports:
	   1) The number of features in B that overlapped the A interval.
	   2) The number of bases in A that had non-zero coverage.
	   3) The length of the entry in A.
	   4) The fraction of bases in A that had non-zero coverage.

## Extracting alignments over specific intervals

You have competing collaborators who want your full dataset to study DPCK and a suite of other genes. You're not interested in sharing the entire dataset before publication, but still want to collaborate. You decide to share with them only the read alignments over the DPCK1 gene. How do you extract a region of aligned reads from a .bam alignment file?

Samtools can do this for you. The DPCK gene is start and stop codons are at 40716757 and 43112247, respectively. Let's isolate reads just for this region and name it DPCK1.sam. By default, samtools will output a SAM file. 

    samtools view -h Brca1Reads_aligned.sorted.dedup.realigned.bam chr17:40716757-43112247 > DPCK1.sam

    head DPCK1.sam

But SAM files are big, we should really compress this into .bam before saving it. What flag in samtools view can you add to output a .bam file instead of the default .sam file?

    samtools view -h

As a relevant aside, now let's package up this data for our collaborators into a tarball, compressed with the bzip2 algorithm.

    mkdir data_for_collaborators
    cp DPCK1.bam chr17.fa data_for_collaborators/
    tar -cjf data_for_collaborators.tar.bz2 data_for_collaborators/

Now we have a nice, small file that we can transfer to our collaborators:

    -rw-rw---- 1 aharkess jlmlab  29M Nov  4 14:26 data_for_collaborators.tar.bz2
    