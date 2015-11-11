In the past 6 labs, you've used a suite of programs and gotten faster at command line tools. Let's combine lessons from a few labs and see how it goes. 

Your collaborators have asked you to identify differentially expressed genes in Zebrafish (_Danio rerio_). They have generated RNAseq for two stages of development, 2 cell embryo and 6 hours post-fertilization. The data are 2x76nt paired-end RNAseq reads derived from polyA-selected mRNA. 

## Write a loop to run FASTQC on each file

First, write a loop to check the quality of the reads using FASTQC. I have given you the data compressed with gzip, which is why they end with .gz. FASTQC can run with gzip-compressed fastq files, so modify your loop accordingly. 

*What you should be asking yourself: Is this data clean? Is there adapter read-through, is there poor sequencing quality?*



