So far in this course, we have been working on one analysis at a time. What if we resequence two, ten, or hundreds of individuals, and want to run the same analysis on all individuals? Since we don't want to individually type commands, then wait, then type again -- we can make the computer work for us while we do something more fun. To do this, we can use loops. Today we will revisit the human brain and adrenal RNAseq from the transcriptomics lab and work on an automation pipeline to save us time. 

The most basic concept of a loop is that we "loop over" multiple things and do something. We will talk about two types of loops today, the "for" loop and the "while" loop. We will implement both types in the bash language in unix. 

## The for loop

Here is a basic skeleton example of a for loop in bash.

    #!/bin/bash
    for i in *.fastq
    do
    SOMETHING INTERESTING
    done
    

Let's start by running FASTQC on the brain and adrenal RNAseq fastq files. These are single end reads (not paired end). Copy the fastqs, chr19.fa, and the gtf gene models into your own directory. I have placed the data here:

    aharkess@zcluster:/home/student/binf4550/data/06.Loops/human_rnaseq$ ls -lh
    total 88M
    -rw-rw---- 1 aharkess jlmlab 7.9M Oct 19 13:56 Adrenal.fq
    -rw-rw---- 1 aharkess jlmlab 6.0M Oct 19 13:56 Brain.fq
    -rw-rw---- 1 aharkess jlmlab  58M Oct 19 13:56 chr19.fa
    -rw-rw---- 1 aharkess jlmlab 5.9M Oct 19 13:57 USCS_hg19_chr19.genes.gtf

To run FASTQC on both of our fastq files in one fell swoop, we can write a loop.

    #!/bin/bash
    for i in *.fq
    do
    echo "Beep boop. I'm a computer running fastqc on file $i"
    /usr/local/fastqc/latest/fastqc $i
    done

## The while loop