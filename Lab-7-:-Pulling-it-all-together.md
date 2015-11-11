In the past 6 labs, you've used a suite of programs and gotten faster at command line tools. Let's combine lessons from a few labs and let you perform your own analysis. For this lab, I will **not** be giving you the full commands to type. Instead, I will provide you a roadmap of programs that you need to execute yourself. You will need to utilize the manual, the help pages and the internet to complete the analysis. 

Your collaborators have asked you to identify differentially expressed genes in Zebrafish (_Danio rerio_). They have generated RNAseq for two stages of development, 2 cell embryo and 6 hours post-fertilization. The data are 2x76nt paired-end RNAseq reads derived from polyA-selected mRNA. These data are NOT strand-specific (fr-unstranded).

![](http://www.zebrafishlab.be/sites/default/files/styles/media_gallery_large/public/embryos-7.jpg)
source: www.zebrafishlab.be

## Get the data and decompress it

I have placed a compressed file with the data on my lab server. **Log onto an interactive node using qlogin** and download this data to your account on the cluster using wget. 

    http://jlmwiki.plantbio.uga.edu/~aharkess/zebrafish_data_from_collaborators.tar.bz2

This file ends in tar.bz2, meaning that it is a tarball (.tar) compressed with bzip2 (.bz2). You have decompressed and unpacked a .tar.bz2 file before -- look in the transcriptome lab to refresh yourself. 

## Write a loop to run FASTQC on each file

First, write a for loop to check the quality of the reads using FASTQC. I have given you the data compressed with gzip, which is why they end with .fastq.gz. FASTQC can run with gzip-compressed fastq files, so modify your loop accordingly. 

**What you should be asking yourself:**

* Is this data clean? 
* Is there adapter read-through? 
* Is there poor sequencing quality?

## Align the reads to the _Danio rerio_ genome using tophat2

Here are steps you need to take to align reads and call differential expression using the Tuxedo suite of programs (specifically, tophat2 and cuffdiff):

1. Index the genome fasta file using bowtie2-build

    /usr/local/bowtie2/latest/bin/bowtie2-build -h

2. Run tophat2, and be sure to include the --GTF flag to specify known transcript locations in your gtf annotation file.

    /usr/local/tophat/latest/bin/tophat2 -h

3. Run cuffdiff to identify differentially expressed genes between the two conditions

    /usr/local/cufflinks/latest/bin/cuffdiff -h
