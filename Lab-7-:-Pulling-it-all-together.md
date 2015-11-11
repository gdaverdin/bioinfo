In the past 6 labs, you've used a suite of programs and gotten faster at command line tools. Let's combine lessons from a few labs and see how it goes. 

Your collaborators have asked you to identify differentially expressed genes in Zebrafish (_Danio rerio_). They have generated RNAseq for two stages of development, 2 cell embryo and 6 hours post-fertilization. The data are 2x76nt paired-end RNAseq reads derived from polyA-selected mRNA. 

![](http://www.zebrafishlab.be/sites/default/files/styles/media_gallery_large/public/embryos-7.jpg)
source: www.zebrafishlab.be

## Get the data and decompress it

I have placed a compressed file with the data on my lab server. Download this data to your account on the cluster using wget. 

    http://jlmwiki.plantbio.uga.edu/~aharkess/zebrafish_data_from_collaborators.tar.bz2

This file ends in tar.bz2, meaning that it is a tarball (.tar) compressed with bzip2 (.bz2). You have decompressed and unpacked a .tar.bz2 file before -- look in the transcriptome lab to refresh yourself. 

## Write a loop to run FASTQC on each file

First, write a loop to check the quality of the reads using FASTQC. I have given you the data compressed with gzip, which is why they end with .gz. FASTQC can run with gzip-compressed fastq files, so modify your loop accordingly. 

**What you should be asking yourself: 

* Is this data clean? 
* Is there adapter read-through? 
* Is there poor sequencing quality?**

## Align the reads to the _Danio rerio_ genome using tophat2




